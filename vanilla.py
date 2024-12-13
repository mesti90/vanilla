#!/usr/bin/python3

#TODO: database??

#Files needed: config.yaml, samples.tsv, nr.dmnd
#Samples.tsv: Sample	Suffix	Reference (e.g. 107-E1	S71	Aci107 for 107-E1_S71)
#config.yaml e.g.:
##threads: 6
##minithreads: 10
##samples: samples.tsv
##fastq: Fastq
##references: Reference
##results: Results
##database: nr




import os
import logging
import time
import yaml
import glob
import csv
import argparse
import multiprocessing as mp
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import pandas as pd
from Bio import GenBank
from Bio.GenBank import Record

def read_args():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-c", help="Config file", default="config.yaml")
	return parser.parse_args()

def init_logging(log_file):
	"""Initialize logging with INFO level."""
	logging.basicConfig(
		format="%(asctime)s. %(levelname)s: %(message)s",
		filename=log_file,
		level=logging.INFO,
	)

def msg(x):
	print(x)
	logging.info(x)

def error(x):
	print(f"Error: {x}")
	logging.error(x)

def mkdir_force(path):
	"""Create a directory if it does not exist."""
	if not os.path.isdir(path):
		os.mkdir(path)

def rm_force(filepath):
	"""Forcefully remove a file if it exists."""
	if os.path.isfile(filepath):
		os.remove(filepath)
		msg(f"Removed file: {filepath}")

def sysexec(command, log_output=True):
	"""Execute a system command and optionally log its output."""
	if log_output:
		msg(f"Executing: {command}")
	try:
		subprocess.run(command, shell=True, check=True)
	except subprocess.CalledProcessError as e:
		error(f"Command failed: {command}\nError: {str(e)}")
		raise

def read_cmd(command):
	"""Execute a shell command and return its output."""
	try:
		result = subprocess.check_output(command, shell=True, universal_newlines=True)
		return result.strip()
	except subprocess.CalledProcessError as e:
		error(f"Error executing command: {command}\n{e}")
		return None

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

def merge_lanes(sample, config, lanes=4):
	"""Merge lane files for a sample."""
	r1_files = sorted(glob.glob(f"{config['fastq']}/{sample}_*_L*_R1_001.fastq.gz"))
	r2_files = sorted(glob.glob(f"{config['fastq']}/{sample}_*_L*_R2_001.fastq.gz"))
	merged_r1 = f"{config['fastq']}/{sample}.R1.fastq.gz"
	merged_r2 = f"{config['fastq']}/{sample}.R2.fastq.gz"

	if len(r1_files) < lanes or len(r2_files) < lanes:
		logging.warning(f"{sample}: Not all lanes downloaded.")
		return

	if os.path.isfile(merged_r1) and os.path.isfile(merged_r2):
		msg(f"{sample}: Merged files already exist.")
		return

	# Merge files
	try:
		sysexec(f"cat {' '.join(r1_files)} > {merged_r1}")
		sysexec(f"cat {' '.join(r2_files)} > {merged_r2}")
		msg(f"{sample}: Merged {len(r1_files)} R1 and {len(r2_files)} R2 files.")
	except Exception as e:
		error(f"Error merging lanes for {sample}: {e}")

def cutadapt(sample, config):
	"""Trim reads using Cutadapt."""
	trimmed_r1 = f"{config['trimmed']}/{sample}.R1.trimmed.fastq.gz"
	trimmed_r2 = f"{config['trimmed']}/{sample}.R2.trimmed.fastq.gz"

	if os.path.isfile(trimmed_r1) and os.path.isfile(trimmed_r2) and not config.get("force"):
		msg(f"{sample}: Trimmed files already exist.")
		return

	try:
		sysexec(
			f"cutadapt --quiet -j {config['threads']} --trim-n --quality-base 33 "
			f"--max-n 0.5 -m 30 --quality-cutoff 30 "
			f"{config['fastq']}/{sample}.R1.fastq.gz "
			f"{config['fastq']}/{sample}.R2.fastq.gz "
			f"-o {trimmed_r1} -p {trimmed_r2}"
		)
		msg(f"{sample}: Trimming completed.")
	except Exception as e:
		error(f"Error trimming reads for {sample}: {e}")


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

def snippy(sample, config):
	"""Run Snippy for variant calling."""
	try:
		tmp_dir = f"{config['snippy']}/{sample.name}_snippy_tmp"
		out_dir = f"{config['snippy']}/{sample.name}_snippy_out"
		mkdir_force(tmp_dir)
		mkdir_force(out_dir)

		sysexec(
			f"snippy --outdir {out_dir} "
			f"--R1 {config['trimmed']}/{sample.name}.R1.trimmed.fastq.gz "
			f"--R2 {config['trimmed']}/{sample.name}.R2.trimmed.fastq.gz "
			f"--cpus {config['threads']} "
			f"--tmpdir {tmp_dir} "
			f"--reference {config['references']}/{sample.reference_name}.gbk --force"
		)
		msg(f"Snippy analysis completed for {sample.name}.")
	except Exception as e:
		error(f"Error running Snippy for {sample.name}: {e}")


class Sample:
	def __init__(self, sample, reference_name, reference_file):
		self.name = sample
		self.reference_name = reference_name
		self.reference_file = reference_file

	def __str__(self):
		return f"Sample: {self.name}, Reference: {self.reference_name}"

	def __repr__(self):
		return f"<Sample(sample={self.name}, reference_name={self.reference_name}, reference_file={self.reference_file})>"

def read_samples(file_path):
	"""Read the sample table and create Sample instances."""
	samples = []
	try:
		# Read the file into a DataFrame, skip comment lines
		df = pd.read_csv(file_path, sep="\t", comment='#')
		
		# Create Sample instances for each row
		for _, row in df.iterrows():
			sample = Sample(
				sample=row['Sample'], 
				reference_name=row['Reference_name'], 
				reference_file=row['Reference_file']
			)
			samples.append(sample)
	except Exception as e:
		print(f"Error reading sample table: {e}")
	
	return samples

def process_sample(sample, config):
	"""Process a single sample through the workflow."""
	try:
		# 1. Trim raw reads
		cutadapt(sample.name, config)

		# 2. Find nucleotide variants
		snippy(sample, config)

		# 3. Identify amino acid variants (custom logic can be added here)
		#TODO
		msg(f"Processing for {sample} completed.")
	except Exception as e:
		error(f"Error processing sample {sample}: {e}")


def predict_orf(base_file, fa, refname, param):
	if os.path.isfile(fa):
		msg(f"{fa} is ready, skipping ORF prediction")
		return
	sysexec(f"getorf -methionine -find {param} -minsize 90 -sequence {base_file} -outseq {fa}")
	convert_orf_ids(fa,refname)

def gunzip_if_necessary(infile, outfile):
	if os.path.isfile(outfile):
		return
	if infile.endswith(".gz"):
		sysexec(f"gunzip -c {infile} > {outfile}")
		return outfile
	return infile

def create_genbank(nt_file, faa_file, gbk_file, reference):
	if os.path.isfile(gbk_file):
		return
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
			genbank_rcd.features.append(ft)
		genbank_records.append(genbank_rcd)
	
	with open(gbk_file,'w') as g:
		for genbank_rcd in genbank_records:
				g.write(str(genbank_rcd))
	msg(f"Ready with {gbk_file}")


def annotate_orfs_with_diamond(fa, diamond_output, config):
	if os.path.isfile(diamond_output):
		return
	sysexec(f"diamond blastp --max-target-seqs 1 --db {config['database']} --threads {config['minithreads']} --out {diamond_output} --outfmt 6 qseqid sseqid qlen slen stitle qstart qend sstart send evalue bitscore length pident qcovhsp scovhsp full_qseq full_sseq --header --query {fa}")

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
	

def predict_and_annotate_orf(reference, config):
	"""Predict ORFs, create GenBank files, and annotate ORFs for a reference genome."""
	refname, reffile = reference
	msg(f"Annotating {refname}")
	wdir = config['references']
	outfaa = f"{wdir}/{refname}.orf.faa"
	outfna = f"{wdir}/{refname}.orf.fna"
	outgbk = f"{wdir}/{refname}.gbk"
	diamond_output = f"{wdir}/{refname}.diamond.top.tsv"
	outbed = f"{wdir}/{refname}.bed"
	base_file = gunzip_if_necessary(reffile, f"{wdir}/{refname}.fna")
	
	for fa, param in [(outfaa, 1,), (outfna, 3,)]:
		predict_orf(base_file, fa, refname, param)
	
	create_genbank(base_file, outfaa, outgbk, refname)
	annotate_orfs_with_diamond(outfaa, diamond_output, config)
	create_annotated_bed(outfna, diamond_output, outbed)
	
	msg(f"Completed annotation for {refname}")




def main():
	# Parse arguments
	args = read_args()

	# Initialize and configure
	config = read_config(args.c)
	log_file = f"{os.path.splitext(__file__)[0]}_{time.strftime('%Y%m%d%H%M%S')}.log"
	init_logging(log_file)

	# Prepare results directory
	mkdir_force(config["results"])
	mkdir_force(config["references"])
	mkdir_force(config["snippy"])

	# Read and process samples
	samples = read_samples(config['samples'])

	# Process references
	references = set((sample.reference_name, sample.reference_file,) for sample in samples)
	msg(f"Processing {len(references)} unique references.")
	for ref in references:
		predict_and_annotate_orf(ref, config)
	#with mp.Pool(config['threads']) as pool:
	#	pool.starmap(predict_and_annotate_orf, [(ref, config) for ref in references])
		
	# Process each sample
	for sample in samples:
		process_sample(sample, config)
		exit()
	#with mp.Pool(config['threads']) as pool:
	#	pool.starmap(process_sample, [(sample, config) for sample in samples])

	# Finalize results
	msg("Pipeline completed successfully.")

if __name__ == "__main__":
	main()
