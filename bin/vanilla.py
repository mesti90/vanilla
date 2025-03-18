#!/usr/bin/python3

#TODO: check if reads are not corrupt; check if trimmed reads are not corrupt

# Standard library imports
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

# Third-party imports
import pandas as pd
import yaml
from Bio import GenBank, SeqIO
from Bio.GenBank import Record
from Bio.Seq import Seq

def read_args():
	"""Parse and return command-line arguments."""
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-c", help="Config file", default="config/vanilla.config.yaml")
	parser.add_argument("-s", help="Samples", default="config/vanilla.samples.tsv")
	return parser.parse_args()

def init_logging(log_file):
	"""Initialize logging with INFO level."""
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

def check_gz_integrity(gz):
	"""Check the integrity of a gzip file."""
	if not os.path.isfile(gz):
		logging.error(f"{gz} does not exist.")
		return False
	res = os.system(f"gzip -t {gz} > /dev/null 2>&1") == 0
	msg(f"{gz} is ok" if res else f"ERROR: {gz} is corrupt")
	return res



def cutadapt(sample):
	"""Trim reads using Cutadapt."""
	trimmed_r1 = f"{config['trimmed']}/{sample}.R1.trimmed.fastq.gz"
	trimmed_r2 = f"{config['trimmed']}/{sample}.R2.trimmed.fastq.gz"
	
	if not config.get("force") and os.path.isfile(trimmed_r1) and os.path.isfile(trimmed_r2):
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

def snippy(sample):
	"""Run Snippy for variant calling."""
	try:
		tmp_dir = os.path.join(config['snippy'],f"{sample.name}_snippy_tmp")
		out_dir = os.path.join(config['snippy'],f"{sample.name}_snippy_out")
		outfile = os.path.join(out_dir, "snps.csv")
		resultfile = os.path.join(config['results'],f"{sample.name}.variants.csv")
		logfile = os.path.join(tmp_dir, "snippy.log")
		errfile = os.path.join(tmp_dir, "snippy.err")
		if os.path.isfile(outfile):
			msg(f"{sample.name}: snippy is skipped, as {outfile} already exists")
		else:
			mkdir_force(tmp_dir)
			mkdir_force(out_dir)

			sysexec(
				f"snippy --outdir {out_dir} "
				f"--R1 {config['trimmed']}/{sample.name}.R1.trimmed.fastq.gz "
				f"--R2 {config['trimmed']}/{sample.name}.R2.trimmed.fastq.gz "
				f"--cpus {config['minithreads']} "
				f"--tmpdir {tmp_dir} "
				f"--reference {config['references']}/{sample.reference_name}.gbk --force > {logfile} 2> {errfile}"
			)
		msg(f"Snippy analysis completed for {sample.name}. Results: {outfile}")
	except Exception as e:
		error(f"Error running Snippy for {sample.name}: {e}")


def annotate_snippy_result(sample):
	snippy_out = os.path.join(config['snippy'],f"{sample.name}_snippy_out", "snps.csv")
	annotated = os.path.join(config["results"], f"{sample.name}.variants.tsv")
	ref_table = self.re

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

def process_sample(sample):
	"""Process a single sample through the workflow."""
	try:
		# 1. Trim raw reads
		cutadapt(sample.name)

		# 2. Find nucleotide and amino acid variants
		snippy(sample)
		
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
		return outfile
	if infile.endswith(".gz"):
		sysexec(f"gunzip -c {infile} > {outfile}")
		return outfile
	return infile

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


def annotate_orfs_with_diamond(fa, diamond_output):
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
	

def predict_and_annotate_orf(reference):
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
	
	annotate_orfs_with_diamond(outfaa, diamond_output)
	create_annotated_bed(outfna, diamond_output, outbed)
	create_genbank(base_file, outfaa, diamond_output, outgbk, refname)
	
	msg(f"Completed annotation for {refname}")

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




def main():
	global config
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
	samples = read_samples(args.s)

	# Predict and annotate ORFs for references
	references = set((sample.reference_name, sample.reference_file,) for sample in samples)
	msg(f"Processing {len(references)} unique references.")
	with mp.Pool(config['threads']) as pool:
		pool.map(predict_and_annotate_orf, references)
	
	# Parallel processing for samples
	with mp.Pool(config['threads']) as pool:
		pool.map(process_sample, samples)

	concat_results(samples, config['variants'])
	# Finalize results
	msg("Pipeline completed successfully.")

if __name__ == "__main__":
	main()
