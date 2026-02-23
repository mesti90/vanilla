#!/usr/bin/python3

"""
This pipeline performs the following steps for each sample:


1. **ORF prediction** using EMBOSS getorf.
2. **Protein annotation** using DIAMOND.
3. **GenBank creation** for predicted ORFs.
"""


#######################
# imports
#######################
import argparse
import os
import subprocess
from collections import defaultdict
import pandas as pd
from Bio import GenBank, SeqIO
from Bio.GenBank import Record
from Bio.Seq import Seq
import shlex
from dataclasses import dataclass
from functools import cached_property
import pdb
pdb.set_trace = lambda *args, **kwargs: None
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

#######################
# singularity utilities
#######################
CONTAINERS = {
	"emboss": "/home/vasarhelyib/containers/staphb-emboss.6.6.0.sif",
	"diamond": "/home/vasarhelyib/containers/buchfink-diamond.2.1.11.sif",
}
_client_bind = "/node8_R10,/node8_data,/node10_R10,/scratch"

#######################
# config and logging
#######################
def msg(x):
	"""Print and log the message"""
	print(x)

def parse_args():
	"""Parse and return command-line arguments."""
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-s", help="Samples", default="samples_for_annotation.tsv")
	parser.add_argument("-n","--threads", help="Max threads", type=int, default=30)
	parser.add_argument("--subthreads", help="Max threads under a main thread", type=int, default=2)
	parser.add_argument("--db", help="Annotation database", default="/node8_R10/vasarhelyib/MinION_results/KBAGDK_20260129/Variant_analysis_20260211/Reference_Staphy/GCF_000013425.1.dmnd")
	args = parser.parse_args()
	args.s = Path(args.s)
	args.db = Path(args.db)
	return args


#######################
# Generic utils
#######################
def singularity_exec(container_name, cmd):
	"""
	Execute a command inside a Singularity container using subprocess.
	
	cmd: str or list, the command to run inside the container
	"""
	container = CONTAINERS[container_name]
	
	if isinstance(cmd, str):
		cmd = cmd.strip()
		cmd_list = ["singularity", "exec", "--bind", _client_bind, container] + cmd.split()
	else:
		cmd_list = ["singularity", "exec", "--bind", _client_bind, container] + cmd

	msg(f"Running: {' '.join(cmd_list)}")
	
	# Run the command and check for errors
	result = subprocess.run(cmd_list, capture_output=True, text=True)
	
	if result.returncode != 0:
		raise RuntimeError(
			f"Command failed in container {container_name}:\n{result.stderr}"
		)
	return result.stdout

#######################
# Sample management
#######################

@dataclass
class Sample:
	name: str
	assembly: Path
	gbk: Path
	
	def __post_init__(self):
		# Ensure the directory for GBK/orf outputs exists
		self.gbk.parent.mkdir(parents=True, exist_ok=True)

	@cached_property
	def orf(self) -> Path:
		"""Return path to ORF protein file (.faa) derived from GenBank filename."""
		return self.gbk.with_suffix(".orf.faa")

	@cached_property
	def orf_fna(self) -> Path:
		"""Return path to ORF nt file (.fna) derived from GenBank filename."""
		return self.gbk.with_suffix(".orf.fna")

	@cached_property
	def diamond(self):
		return self.gbk.with_suffix(".diamond.tsv")


	@cached_property
	def bed(self) -> Path:
		"""Return path to BED file derived from GenBank filename."""
		return self.gbk.with_suffix(".bed")

	def __str__(self):
		return f"Sample: {self.name}"
	
	def __repr__(self):
		return f"<{self.__str__}>"


def read_sample_table(file_path):
	"""Read the sample table and create Sample instances."""
	samples = []
	try:
		# Read the file into a DataFrame, skip comment lines
		df = pd.read_csv(file_path, sep="\t", comment='#').fillna("")
		required_columns = set(['Sample','Assembly','Genbank'])
		present = set(df.columns)
		if missing := required_columns - present:
			raise ValueError(f"Missing columns: {'; '.join(missing)}")
		# Create Sample instances for each row
		for _, row in df.iterrows():
			if not row['Sample'].startswith("#"):
				sample = Sample(name=row['Sample'], assembly=Path(row['Assembly']), gbk=Path(row['Genbank']))
				samples.append(sample)
	except Exception as e:
		print(f"Error reading sample table: {e}")
	
	return samples



#######################
# ORF prediction and annotation
#######################

def predict_orf(sample, param, outfile):
	if outfile.exists():
		msg(f"{outfile} is ready, skipping ORF prediction")
		return
	singularity_exec("emboss", f"getorf -methionine -find {param} -minsize 90 -sequence {sample.assembly} -outseq {outfile}")
	convert_orf_ids(sample,outfile)


def convert_orf_ids(sample,orf_file):
	"""Standardize ORF IDs in a FASTA file."""
	updated_records = []
	for rcd in SeqIO.parse(orf_file,"fasta"):
		block = rcd.description.strip().split()
		ctg, orfnum = block[0].rsplit("_",1)
		start,end = sorted([int(block[1].strip("[")),int(block[3].strip("]"))])
		ori = "-" if rcd.description.strip().endswith("(REVERSE SENSE)") else "+"
		orfid = f"{sample.name}|{ctg}|ORF_{orfnum}|{ori}|{start}-{end}"
		updated_records.append(">{}\n{}\n".format(orfid,rcd.seq))
		
	with open(orf_file, "w") as g:
		g.writelines(updated_records)

def annotate_orfs_with_diamond(sample, args):
	if sample.diamond.exists():
		return
	cmd = f"diamond blastp --max-target-seqs 1 --db {args.db} --threads {args.subthreads} --out {sample.diamond} --outfmt 6 qseqid sseqid qlen slen stitle qstart qend sstart send evalue bitscore length pident qcovhsp scovhsp full_qseq full_sseq --header --query {sample.orf}"
	singularity_exec("diamond", cmd)


def create_genbank(sample):
	if sample.gbk.exists():
		return
	df = pd.read_csv(sample.diamond, sep="\t", comment='#', header=None)
	products = {row[0]: row[4] for _,row in df.iterrows()}
	msg(f"Reading {sample.orf_fna} and {sample.orf}")
	fa_nt = SeqIO.index(sample.assembly,"fasta")
	fa_aa = SeqIO.parse(sample.orf,"fasta")
	contigs = list(fa_nt.keys())
	orfs = {ctg:[] for ctg in contigs}
	
	pdb.set_trace()
	
	for rcd in fa_aa:
		ctg = rcd.id.split("|")[1]
		orfs[ctg].append([rcd.id,rcd.seq])
	
	msg(f"Creating genbank file for {sample.name}")

	genbank_records = []
	for ctg in contigs:
		genbank_rcd = Record.Record()
		genbank_rcd.sequence = str(fa_nt[ctg].seq)
		genbank_rcd.locus = f"{sample.name}|{ctg}"
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
	
	with open(sample.gbk,'w') as g:
		for genbank_rcd in genbank_records:
				g.write(f"{genbank_rcd}\n")
	msg(f"Ready with {sample.gbk}")




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
	

def create_annotated_bed(sample):
	bed, orfid_list = get_bed_entries(sample.orf_fna)
	annotation = defaultdict(lambda:"-")
	with open(sample.diamond) as f:
		for line in f:
			line = line.strip()
			if line and not line.startswith("#"):
				block = line.split("\t")
				annotation[block[0]]= block[4]
	with open(sample.bed,"w") as g:
		for orfid in orfid_list:
			g.write("\t".join(bed[orfid]) +"\t{}\n".format(annotation[orfid]))
	

def annot_sample(arglist):
	"""Predict ORFs, create GenBank files, and annotate ORFs"""
	sample, args = arglist
	if sample.gbk.exists():
		return
	msg(f"Annotating {sample.name}")
	
	predict_orf(sample, param=1, outfile=sample.orf)
	predict_orf(sample, param=3, outfile=sample.orf_fna)
	annotate_orfs_with_diamond(sample, args)
	create_annotated_bed(sample)
	create_genbank(sample)
	
	msg(f"Completed annotation for {sample.name}")

def main():
	# Parse arguments
	args = parse_args()

	samples = read_sample_table(args.s)
	
	msg(f"Processing {len(samples)} unique samples.")
	#for sample in samples:
	#	annot_sample((sample, args,))
	with ProcessPoolExecutor(max_workers=args.threads) as executor:
		list(executor.map(annot_sample, [(sample, args,) for sample in samples]))
	
	msg("Pipeline completed successfully.")

if __name__ == "__main__":
	main()
