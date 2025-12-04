#!/usr/bin/python3


"""
Comprehensive Bioinformatics Pipeline
------------------------------------
This pipeline performs the following steps for each sample:


1. **Read trimming** using Cutadapt.
2. **Variant calling** using Snippy.
3. **ORF prediction** using EMBOSS getorf.
4. **Protein annotation** using DIAMOND.
5. **GenBank creation** for predicted ORFs.
6. **Aggregation of results** across samples.


It also processes each *reference genome* used across samples by predicting ORFs,
annotating proteins, generating GenBank files, and preparing BED tracks.


Key components:
- YAML configuration file specifying paths and parameters.
- Samples TSV listing sample name, reference identifier, and reference file.
- External tools required: Cutadapt, Snippy, GetOrf, Diamond, gzip.


The code is structured for readability and modularity. Each major action is handled
by one or more helper functions, with extensive logging included.
"""


#######################
# imports
#######################
import argparse
import csv
import glob
import logging
import multiprocessing as mp
import os
import subprocess
import time
from collections import defaultdict
import shutil
import pandas as pd
import yaml
from Bio import GenBank, SeqIO
from Bio.GenBank import Record
from Bio.Seq import Seq
from spython.main import Client
import shlex
from dataclasses import dataclass, field
import pdb
pdb.set_trace = lambda *args, **kwargs: None
from concurrent.futures import ProcessPoolExecutor


#######################
# spython utilities
#######################
CONTAINERS = {
	"emboss": "/home/vasarhelyib/containers/staphb-emboss.6.6.0.sif",
	"diamond": "/home/vasarhelyib/containers/buchfink-diamond.2.1.11.sif",
	"cutadapt": "/home/vasarhelyib/containers/mesti90-cutadapt.5.1.sif",
	"snippy": "/home/vasarhelyib/containers/staphb-snippy-4.6.0-SC2.img",
}
_client_options = ["--bind", "/node8_R10,/node8_data,/node10_R10,/scratch"]
INSTANCES = {key: Client.instance(value, options=_client_options) for key, value in CONTAINERS.items()}

lock = mp.Lock()

#######################
# config and logging
#######################
# Logging infrastructure
def init_logging(log_file):
	"""Initialize logging with INFO level."""
	os.makedirs(os.path.dirname(log_file), exist_ok=True)
	logging.basicConfig(
		format="%(asctime)s. %(levelname)s: %(message)s",
		filename=log_file,
		level=logging.INFO,
	)

def msg(x):
	"""Print and log the message"""
	print(x)
	logging.info(x)

def error(x):
	"""Prints an error"""
	print(f"Error: {x}")
	logging.error(x)


def read_args():
	"""Parse and return command-line arguments."""
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-c", help="Config file", default="config/vanilla.config.yaml")
	parser.add_argument("-s", help="Samples", default="config/vanilla.samples.tsv")
	parser.add_argument("-a", "--make-annotation", action="store_true", help="Generate ORF prediction + annotation if missing Reference_gbk")
	parser.add_argument("-f", "--force", action="store_true", help="Overwrite existing files")
	return parser.parse_args()


def read_config(config_file):
	"""Load configuration from a YAML file."""
	try:
		with open(config_file) as f:
			config = yaml.safe_load(f)
		# Default values
		config.setdefault("project", "output")
		config.setdefault("results", "Results")
		for key, value in config.items():
			if isinstance(value, str):
				config[key] = value.format(project=config["project"], results=config["results"])
		return config
	except Exception as e:
		error(f"Error reading config file {config_file}: {e}")
		raise



#######################
# Generic utils
#######################
def spython_exec(container_name, cmd):
	"""Execute a command inside a Singularity container instance."""
	cmdlist = cmd if isinstance(cmd, list) else shlex.split(str(cmd))
	return Client.execute(INSTANCES[container_name], cmdlist)

def spython_run(container_name, cmd):
	"""Execute a command inside a Singularity container instance."""
	cmdlist = cmd if isinstance(cmd, list) else shlex.split(str(cmd))
	return Client.run(INSTANCES[container_name], cmdlist)


def mkdir_force(path):
	"""Create a directory if it does not exist."""
	os.makedirs(path, exist_ok=True)

def rm_force(filepath):
	"""Forcefully remove a file if it exists."""
	if os.path.isfile(filepath):
		os.remove(filepath)
		msg(f"Removed file: {filepath}")

def check_gz_integrity(gz):
	"""Check the integrity of a gzip file."""
	if not os.path.isfile(gz):
		logging.error(f"{gz} does not exist.")
		return False
	res = os.system(f"gzip -t {gz} > /dev/null 2>&1") == 0
	msg(f"{gz} is ok" if res else f"ERROR: {gz} is corrupt")
	return res

def gunzip_if_necessary(infile, outfile):
	if os.path.isfile(outfile):
		return outfile
	if infile.endswith(".gz"):
		subprocess.run(f"gunzip -c {infile} > {outfile}", shell=True, check=True)
		return outfile
	return infile

#######################
# Sample management
#######################

@dataclass
class Sample:
	name: str
	r1: str
	r2: str
	fa: str
	reference_name: str
	reference_fasta: str
	reference_gbk: str
	config: dict
	# Derived fields
	trimmed1: str = field(init=False)
	trimmed2: str = field(init=False)

	def __post_init__(self):
		"""Determine sample type and output trimmed file paths."""
		# contig-based sample
		pdb.set_trace()
		if self.fa:
			return
		if self.r1 and self.r2:
			self.trimmed1 = f"{self.config['trimmed']}/{self.name}.R1.trimmed.fastq.gz"
			self.trimmed2 = f"{self.config['trimmed']}/{self.name}.R2.trimmed.fastq.gz"
			return
		# missing input
		raise ValueError(f"Invalid sample '{self.name}': must provide contigs or both R1 and R2 reads.")

	def __str__(self):
		return f"Sample: {self.name}, Reference: {self.reference_name}"
	
	def __repr__(self):
		return f"<Sample(name={self.name}, reference={self.reference_name})>"


def read_samples(file_path,config):
	"""Read the sample table and create Sample instances."""
	samples = []
	try:
		# Read the file into a DataFrame, skip comment lines
		df = pd.read_csv(file_path, sep="\t", comment='#').fillna("")
		required_columns = set(['Sample','Reference_name'])
		present = set(df.columns)
		if missing := required_columns - present:
			raise ValueError(f"Missing columns: '; '.join({missing})")
		# Create Sample instances for each row
		for _, row in df.iterrows():
			sample = Sample(
				name=row['Sample'],
				r1=row.get('R1', ""),
				r2=row.get('R2', ""),
				fa=row.get('contigs', ""),
				reference_name=row['Reference_name'],
				reference_fasta=row.get('Reference_fasta', ""),
				reference_gbk=row.get('Reference_gbk', ""),
				config=config
			)

			samples.append(sample)
	except Exception as e:
		print(f"Error reading sample table: {e}")
	
	return samples


#######################
# Main functions for Workflow
#######################


#######################
# Trimming
#######################
def cutadapt(sample):
	"""Trim reads using Cutadapt."""
	
	if not config.get("force") and os.path.isfile(sample.trimmed1) and os.path.isfile(sample.trimmed2):
		msg(f"{sample}: Trimmed files already exist.")
		return
	
	with lock:
		parent_dir = os.path.dirname(sample.trimmed1)
		# Create parent directory if it doesn't exist
		os.makedirs(parent_dir, exist_ok=True)
	
	try:
		spython_run('cutadapt', f"cutadapt -j {config['threads']} --trim-n --quality-base 33 --max-n 0.5 -m 30 --quality-cutoff 30 {sample.r1} {sample.r2} -o {sample.trimmed1} -p {sample.trimmed2}")
		msg(f"{sample}: Trimming completed.")
	except Exception as e:
		error(f"Error trimming reads for {sample}: {e}")

#######################
# ORF prediction and annotation
#######################

def predict_orf(base_file, fa, refname, param):
	if os.path.isfile(fa):
		msg(f"{fa} is ready, skipping ORF prediction")
		return
	spython_exec("emboss", f"getorf -methionine -find {param} -minsize 90 -sequence {base_file} -outseq {fa}")
	convert_orf_ids(fa,refname)


def convert_orf_ids(fname,reference):
	"""Standardize ORF IDs in a FASTA file."""
	updated_records = []
	for rcd in SeqIO.parse(fname,"fasta"):
		block = rcd.description.strip().split()
		ctg, orfnum = block[0].rsplit("_",1)
		start,end = sorted([int(block[1].strip("[")),int(block[3].strip("]"))])
		ori = "-" if rcd.description.strip().endswith("(REVERSE SENSE)") else "+"
		orfid = f"{reference}|{ctg}|ORF_{orfnum}|{ori}|{start}-{end}"
		updated_records.append(">{}\n{}\n".format(orfid,rcd.seq))
		
	with open(fname, "w") as g:
		g.writelines(updated_records)

def annotate_orfs_with_diamond(fa, diamond_output):
	if os.path.isfile(diamond_output):
		return
	cmd = f"diamond blastp --max-target-seqs 1 --db {config['database']} --threads {config['minithreads']} --out {diamond_output} --outfmt 6 qseqid sseqid qlen slen stitle qstart qend sstart send evalue bitscore length pident qcovhsp scovhsp full_qseq full_sseq --header --query {fa}"
	msg(cmd)
	spython_exec("diamond", cmd)


def create_genbank(nt_file, faa_file, diamond_output, gbk_file, reference):
	if os.path.isfile(gbk_file):
		return
	df = pd.read_csv(diamond_output, sep="\t", comment='#', header=None)
	products = {row[0]: row[4] for _,row in df.iterrows()}
	msg(f"Reading {nt_file} and {faa_file}")
	fa_nt = SeqIO.index(nt_file,"fasta")
	fa_aa = SeqIO.parse(faa_file,"fasta")
	contigs = list(fa_nt.keys())
	orfs = {ctg:[] for ctg in contigs}
	
	for rcd in fa_aa:
		ctg = rcd.id.split("|")[1]
		orfs[ctg].append([rcd.id,rcd.seq])
	
	msg(f"Creating genbank file for {reference}")

	genbank_records = []
	for ctg in contigs:
		genbank_rcd = Record.Record()
		genbank_rcd.sequence = str(fa_nt[ctg].seq)
		genbank_rcd.locus = f"{reference}|{ctg}"
		genbank_rcd.size = len(fa_nt[ctg].seq)
		genbank_rcd.residue_type = "DNA"
		for orf in orfs[ctg]:
			block = orf[0].split("|")
			start,end = block[4].split("-")
			location = f"{start}..{end}" if block[3] == "+" else f"complement({start}..{end})"
			ft = Record.Feature("CDS",location)
			ft.qualifiers.append(Record.Qualifier("/gene=",f'"{orf[0]}"'))
			ft.qualifiers.append(Record.Qualifier("/translation=",f'"{orf[1]}"'))
			if orf[0] in products:
				ft.qualifiers.append(Record.Qualifier("/product=",f'"{products[orf[0]]}"'))
			genbank_rcd.features.append(ft)
		genbank_records.append(genbank_rcd)
	
	with open(gbk_file,'w') as g:
		for genbank_rcd in genbank_records:
				g.write(str(genbank_rcd))
	msg(f"Ready with {gbk_file}")




def get_bed_entries(fname):
	"""Generate BED entries and a list of ORF IDs from a FASTA file."""
	bed = {}
	orfid_list = []
	for rcd in SeqIO.parse(fname, "fasta"):
		_, ctg, _, ori, coords = rcd.id.split("|")
		start, end = map(int, coords.split("-"))
		bed[rcd.id] = [ctg, str(start - 1), str(end), rcd.id, ori]
		orfid_list.append(rcd.id)
	return bed, orfid_list
	

def create_annotated_bed(orf_fna, diamond_output, outbed):
	bed, orfid_list = get_bed_entries(orf_fna)
	annotation = defaultdict(lambda:"-")
	with open(diamond_output) as f:
		for line in f:
			line = line.strip()
			if line and not line.startswith("#"):
				block = line.split("\t")
				annotation[block[0]]= block[4]
	with open(outbed,"w") as g:
		for orfid in orfid_list:
			g.write("\t".join(bed[orfid]) +"\t{}\n".format(annotation[orfid]))
	

def predict_and_annotate_orf(refname, reffile, outgbk):
	"""Predict ORFs, create GenBank files, and annotate ORFs for a reference genome."""
	msg(f"Annotating {refname}")
	wdir = config['references']
	outfaa = f"{wdir}/{refname}.orf.faa"
	outfna = f"{wdir}/{refname}.orf.fna"
	#outgbk = f"{wdir}/{refname}.gbk"
	diamond_output = f"{wdir}/{refname}.diamond.top.tsv"
	outbed = f"{wdir}/{refname}.bed"
	base_file = gunzip_if_necessary(reffile, f"{wdir}/{refname}.fna")
	
	for fa, param in [(outfaa, 1,), (outfna, 3,)]:
		predict_orf(base_file, fa, refname, param)
	
	annotate_orfs_with_diamond(outfaa, diamond_output)
	create_annotated_bed(outfna, diamond_output, outbed)
	create_genbank(base_file, outfaa, diamond_output, outgbk, refname)
	
	msg(f"Completed annotation for {refname}")


#######################
# Variant calling
#######################
def snippy(sample):
	"""Run Snippy for variant calling."""
	try:
		tmp_dir = os.path.join(config['snippy'],f"{sample.name}_snippy_tmp")
		out_dir = os.path.join(config['snippy'],f"{sample.name}_snippy_out")
		outfile = os.path.join(out_dir, "snps.csv")
		resultfile = os.path.join(config['results'],f"{sample.name}.variants.csv")
		logfile = os.path.join(tmp_dir, "snippy.log")
		errfile = os.path.join(tmp_dir, "snippy.err")
		if not config['force'] and os.path.isfile(outfile):
			msg(f"{sample.name}: snippy is skipped, as {outfile} already exists")
		else:
			mkdir_force(tmp_dir)
			mkdir_force(out_dir)
			
			reference_input = sample.reference_gbk if sample.reference_gbk else sample.reference_fasta
			snippy_cmd = f"snippy --outdir {out_dir} --cpus {config['minithreads']} --tmpdir {tmp_dir} --reference {reference_input} --force --quiet "
			if sample.fa:
				source_txt = f" --ctgs {sample.fa} "
			elif sample.trimmed1 and sample.trimmed2:
				source_txt = f" --R1 {sample.trimmed1} --R2 {sample.trimmed2} "
			else:
				msg(f"Error with {sample.name}: no reads and no contigs were seen")
				raise Exception
			snippy_cmd += source_txt
			msg(snippy_cmd)
			spython_run("snippy", snippy_cmd)
		msg(f"Snippy analysis completed for {sample.name}. Results: {outfile}")
	except Exception as e:
		error(f"Error running Snippy for {sample.name}: {e}")



#######################
# Process samples
#######################
def process_sample(sample):
	"""Process a single sample through the workflow."""
	try:
		# 1. Trim raw reads
		if not sample.fa:
			cutadapt(sample)

		# 2. Find nucleotide and amino acid variants
		snippy(sample)
		
		msg(f"Processing for {sample} completed.")
	except Exception as e:
		error(f"Error processing sample {sample}: {e}")


#######################
# Concatenating results
#######################
def concat_results(samples, concatenated_file):
	"""Concatenate results from all samples into a single tab-separated file with a unified header.

	Args:
		samples (list): List of Sample objects.
		concatenated_file (str): Path to the output concatenated file.
	"""
	header_written = False  # Flag to track if the header has been written

	with open(concatenated_file, "w", newline='') as outfile:
		tsv_writer = csv.writer(outfile, delimiter='\t')

		for sample in samples:
			out_dir = os.path.join(config['snippy'], f"{sample.name}_snippy_out")
			outfile_path = os.path.join(out_dir, "snps.csv")

			if os.path.isfile(outfile_path):
				with open(outfile_path, "r", newline='') as infile:
					csv_reader = csv.reader(infile)
					for i, row in enumerate(csv_reader):
						if i == 0:  # Handle the header row
							if not header_written:
								tsv_writer.writerow(["Sample"] + row)  # Write the header with "Sample" column
								header_written = True
						else:
							tsv_writer.writerow([sample.name] + row)  # Write the data row with "Sample" column
			else:
				msg(f"Warning: File {outfile_path} does not exist for sample {sample.name}.")
	
	msg(f"Concatenated results saved to {concatenated_file}.")


def process_reference(arglist):
	reference_name, fasta, gbk, make_annotation = arglist
	if os.path.isfile(gbk):
		msg(f"{reference_name}: Using existing annotation {gbk}")
	elif make_annotation:
		msg(f"{reference_name}: Creating annotation (gbk missing)")
		predict_and_annotate_orf(reference_name, fasta, gbk)
	else:
		msg(f"{reference_name}: No gbk annotation and --make-annotation not set; "
			f"snippy will use FASTA only as reference")

def main():
	global config
	# Parse arguments
	args = read_args()

	# Initialize and configure
	config = read_config(args.c)
	log_file = f"Log/{os.path.splitext(os.path.basename(__file__))[0]}_{time.strftime('%Y%m%d%H%M%S')}.log"
	init_logging(log_file)
	
	config['force'] = args.force
	
	# Prepare results directory
	mkdir_force(config["results"])
	mkdir_force(config["references"])
	mkdir_force(config["snippy"])

	# Read and process samples
	samples = read_samples(args.s, config)

	# Predict and annotate ORFs for references
	references = set((sample.reference_name, sample.reference_fasta, sample.reference_gbk) for sample in samples)
	msg(f"Processing {len(references)} unique references.")
	with ProcessPoolExecutor(max_workers=4) as executor:
		for result in executor.map(process_reference, [(reference_name, fasta, gbk, args.make_annotation,) for  reference_name, fasta, gbk in references]):
			msg(result)
	
	# Parallel processing for samples
	with mp.Pool(config['threads']) as pool:
		pool.map(process_sample, samples)

	concat_results(samples, config['variants'])
	# Finalize results
	msg("Pipeline completed successfully.")

if __name__ == "__main__":
	main()
