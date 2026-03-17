#!/usr/bin/env python3
import pdb
import subprocess
import shlex
from datetime import datetime
from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import random
import threading
import pandas as pd
from Bio import SeqIO
import gzip

lock = threading.Lock()

############################################
# CONTAINERS
############################################

KEY_COLS = ['CHROM', 'POS', 'TYPE', 'REF', 'ALT']

class Containers:
	"""Holds paths to all Singularity container images."""
	def __init__(self, base_dir: Path = Path("/node8_R10/vasarhelyib/containers")):
		self.FILT_LONG_IMG = base_dir / "staphb-filtlong.0.3.1.sif"
		self.UNICYCLER_IMG = base_dir / "staphb-unicycler.0.5.1.sif"
		self.RAVEN_IMG = base_dir / "staphb-raven.1.8.3.sif"
		self.FLYE_IMG = base_dir / "staphb-flye.2.9.6.sif"
		self.SNIPPY_IMG = base_dir / "staphb-snippy-4.6.0-SC2.img"

	@property
	def all_images(self):
		return [v for k, v in vars(self).items() if k.endswith("_IMG")]

class Tool:
	def __init__(self, name: str, container_attr: str, run_function):
		"""
		name: Human-readable name
		container_attr: Attribute name in Containers (e.g. "UNICYCLER_IMG")
		run_function: Function to execute (e.g. run_unicycler)
		"""
		self.name = name
		self.container_attr = container_attr
		self.run_function = run_function

	def run(self, sample, args, containers, *extra):
		container = getattr(containers, self.container_attr)
		self.run_function(sample, *extra, args, container)

############################################
# ARGUMENTS
############################################

def parse_args():
	parser = argparse.ArgumentParser(description="Pipeline for genome assembly and variant analysis using Filtlong, Unicycler, Raven, Flye, and Snippy.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them.")
	parser.add_argument("-s", "--sample_table", default="samples_for_vanillong.tsv",
						help="TSV file with sample info: Sample, Reference_gbk, Raw_reads.")
	parser.add_argument("-n", "--threads", type=int, default=20, help="Number of threads per assembly tool.")
	parser.add_argument("-wd", "--workdir", default="vanillong_work", help="Working directory for outputs.")
	parser.add_argument("--subthreads", type=int, default=3, help="Maximum number of samples to process concurrently.")
	parser.add_argument("--binds", default="/node8_R10,/node10_R10,/scratch", help="These folders will be binded for singularity")
	parser.add_argument("--stats", default="variant_statistics.tsv",help="The variant statistics will be printed to this file")
	parser.add_argument("--all_variants", default="all_variants.tsv", help="The variants will be printed here")
	
	parser.add_argument("--filtlong_max_bases", default="500mb", help="This arg will be passed to filtlong")
	
	parser.add_argument("--test", action="store_true")
	parser.add_argument("--call_variants", action="store_true", help="Calling variants")
	parser.add_argument("--create_stats", action="store_true", help="Compute stats from existing snps.csv files.")
	parser.add_argument("--combine_variants", action="store_true", help="Combine per-sample variant tables into one file.")
	args = parser.parse_args()
	args.workdir = Path(args.workdir.strip())
	args.sample_table = Path(args.sample_table.strip())
	args.stats = args.workdir / args.stats
	args.all_variants = args.workdir / args.all_variants
	args.errfile = args.workdir / "vanillong.err"
	return args

############################################
# SAMPLE
############################################


class Sample:
	"""Holds sample-specific data, output directory, and all intermediate results."""
	def __init__(self, name: str, ref_gbk: Path, raw_reads: Path, workdir: Path):
		self.name = name.strip()
		self.ref_gbk = ref_gbk
		self.raw_reads = raw_reads
		self.sample_dir = workdir / self.name

		# Assembly / intermediate files
		self.filtlong_fastq: Path = self.sample_dir / f"{self.name}.filtlong.fastq"
		self.unicycler_dir : Path = self.sample_dir / "unicycler"
		self.unicycler_assembly : Path = self.unicycler_dir / "assembly.fasta"
		self.raven_assembly : Path = self.sample_dir / f"{self.name}.raven.fna"
		self.flye_dir : Path = self.sample_dir / "flye"
		self.flye_assembly : Path =  self.flye_dir / "assembly.fasta"

		# Snippy output dirs
		self.snippy_unicycler_dir : Path = self.sample_dir / "snippy_unicycler"
		self.snippy_raven_dir : Path = self.sample_dir / "snippy_raven"
		self.snippy_flye_dir : Path = self.sample_dir / "snippy_flye"
		self.snippy_raw_dir: Path = self.sample_dir / "snippy_raw"
		
		#Variants
		self.variants = workdir / f"{self.name}.variants.tsv"

	def summary(self):
		"""Print a summary of all intermediate files and results."""
		print(f"Sample: {self.name}")
		print(f"  Filtlong: {self.filtlong_fastq}")
		print(f"  Unicycler assembly: {self.unicycler_assembly}")
		print(f"  Raven assembly: {self.raven_assembly}")
		print(f"  Flye assembly: {self.flye_assembly}")
		print(f"  Snippy Unicycler: {self.snippy_unicycler_dir}")
		print(f"  Snippy Raven: {self.snippy_raven_dir}")
		print(f"  Snippy Flye: {self.snippy_flye_dir}")
		print(f"  Snippy Raw reads: {self.snippy_raw_dir}")
		print("")
############################################
# UTILITIES
############################################

COLOR_YELLOW = "\033[93m"
COLOR_RED = "\033[91m"
COLOR_RESET = "\033[0m"


def print_prefix(prefix, color=None):
	"""Print a message prefix with optional color."""
	ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	pre = f"[{ts}] [{prefix}]"
	if color:
		pre = f"{color}{pre}{COLOR_RESET}"
	return pre

def print_skip(msg):
	"""Print a timestamped SKIP message in yellow."""
	print(f"{print_prefix('SKIP', COLOR_YELLOW)} {msg}")

def print_error(msg):
	"""Print a timestamped ERROR message in red."""
	print(f"{print_prefix('ERROR', COLOR_RED)} {msg}")

def print_info(msg):
	print(f"{print_prefix('INFO')} {msg}")

def log_error(message: str, errfile: Path):
	"""Append error messages with timestamp to an error file."""
	ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	with lock:
		with errfile.open("a") as f:
			f.write(f"[{ts}] {message}\n")


def run_cmd(cmd, dry_run: bool, errfile: Path = None, capture_output=False):
	"""Run a command, optionally suppress stdout/stderr and log errors."""
	print(print_prefix("DRY_RUN" if dry_run else "RUN") + " ".join(shlex.quote(c) for c in cmd))
	if dry_run:
		return None

	try:
		if capture_output:
			result = subprocess.run(cmd, check=True, capture_output=True, text=True)
		else:
			result = subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
			if result.stderr and errfile:
				log_error(result.stderr.strip(), errfile)
		return result
	except subprocess.CalledProcessError as e:
		if errfile:
			log_error(f"Command failed: {' '.join(cmd)}\n{e.stderr or e}", errfile)
		raise

def run_cmd_redirect(outfile: Path, cmd, dry_run: bool, errfile: Path = None):
	print(print_prefix("DRY_RUN" if dry_run else "RUN") + " ".join(shlex.quote(c) for c in cmd) + f" > {shlex.quote(str(outfile))}")
	if dry_run:
		return

	try:
		with outfile.open("w") as f:
			subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
	except subprocess.CalledProcessError as e:
		if errfile:
			log_error(f"Command failed: {' '.join(cmd)}\n{e.stderr or e}", errfile)
		raise


def singularity_cmd(img: Path, binds: str, *args):
	return ["singularity", "run", "-B", binds, str(img), *args]

############################################
# TOOL FUNCTIONS
############################################

def non_empty(fname):
	return fname.exists() and fname.stat().st_size > 0

def run_filtlong(sample: Sample, args, container: Path):
	if non_empty(sample.filtlong_fastq):
		print_skip(f"filtlong for {sample.name}, output exists: {sample.filtlong_fastq}")
		return
		
	run_args = [
		"filtlong",
		"--min_length", "1kb",
		"--min_mean_q", "85",
		"--min_window_q", "60",
		"--keep_percent", "80",
		"--target_bases", args.filtlong_max_bases,
		str(sample.raw_reads)
	]
	cmd = singularity_cmd(container, args.binds, *run_args)
	run_cmd_redirect(sample.filtlong_fastq, cmd, args.dry_run)



def run_unicycler(sample: Sample, args, container: Path):
	if non_empty(sample.unicycler_assembly):
		print_skip( f"unicycler for {sample.name}, output exists: {sample.unicycler_assembly}")
		return
		
	sample.unicycler_dir.mkdir(parents=True, exist_ok=True)
	run_args = [
		"unicycler",
		"--long", str(sample.filtlong_fastq),
		"--out", str(sample.unicycler_dir),
		"--keep", "0",
		"--threads", str(args.threads),
		"--mode", "conservative",
		"--verbosity", "0",
	]
	cmd = singularity_cmd(container, args.binds, *run_args)
	run_cmd(cmd, args.dry_run)
	
		

def run_raven(sample: Sample, args, container: Path) :
	if non_empty(sample.raven_assembly):
		print_skip(f"raven for {sample.name}, output exists: {sample.raven_assembly}")
		return
		
	run_args = ["raven", str(sample.filtlong_fastq), "-t", str(args.subthreads),"--disable-checkpoints"]
	cmd = singularity_cmd(container, args.binds, *run_args)
	run_cmd_redirect(sample.raven_assembly, cmd, args.dry_run)
	
def run_flye(sample: Sample, args, container: Path) :
	if non_empty(sample.flye_assembly):
		print_skip(f"flye for {sample.name}, output exists: {sample.flye_assembly}")
		return
		
	sample.flye_dir.mkdir(parents=True, exist_ok=True)
	run_args = ["flye", "--nano-raw", str(sample.filtlong_fastq), "--out-dir", str(sample.flye_dir), "-t", str(args.subthreads)]
	cmd = singularity_cmd(container, args.binds, *run_args)
	try:
		run_cmd(cmd, args.dry_run)
	except Exception as e:
		print_error(f"running flye for {sample.name}: {e}")
		log_error(f"Error running flye for {sample.name}: {e}\nCommand: {' '.join(cmd)}\n\n", args.errfile)

def run_snippy(sample, input_file: Path,  outdir: Path, args, container: Path, mode: str = "ctg"):
	expected_file = outdir / "snps.csv"  # or pick any file that indicates Snippy finished
	if non_empty(expected_file):
		print_skip( f"snippy ({mode}) for {sample.name}, output exists: {expected_file}")
		return
		
	outdir.mkdir(parents=True,exist_ok=True)
	if not args.dry_run and not input_file.exists():
		print(f"{sample.name}: {input_file} is missing, skipping snippy")
		return
	if mode not in ("ctg", "se"):
		print(f"{input_file}: Invalid mode: {mode}. Must be 'ctg' or 'se'.")
		return
	run_args = ["snippy", "--force", "--outdir", str(outdir), "--ref", str(sample.ref_gbk), "--cpus", str(args.subthreads), f"--{mode}", str(input_file)]
	cmd = singularity_cmd(container, args.binds, *run_args)
	try:
		run_cmd(cmd, args.dry_run)
	except Exception as e:
		print_error(f"running snippy ({mode}) for sample {sample.name}, {sample.ref_gbk} on {input_file.name}: {e}")
		log_error(f"Error running snippy ({mode}) for sample {sample.name}, {sample.ref_gbk} on {input_file.name}: {e}\nCommand: {' '.join(cmd)}\n\n", args.errfile)

############################################
# SAMPLE PROCESSING
############################################

def get_stats_and_variants_for_sample(sample: Sample, args, containers: Containers):
	#pdb.set_trace()
	if not sample.raw_reads.is_file():
		raise FileNotFoundError(f"Missing reads: {sample.raw_reads}")
	if not sample.ref_gbk.is_file():
		raise FileNotFoundError(f"Missing reference: {sample.ref_gbk}")

	sample.sample_dir.mkdir(parents=True, exist_ok=True)

	# Step 1: Filtlong
	filtlong_tool = Tool("filtlong", "FILT_LONG_IMG", run_filtlong)
	filtlong_tool.run(sample, args, containers)
	if not args.dry_run and not non_empty(sample.filtlong_fastq):
		print(f"Error: filtlong for {sample.name} failed")
		return

	# Step 2: Assemblies
	assembly_tools = [
		Tool("unicycler", "UNICYCLER_IMG", run_unicycler),
		Tool("raven", "RAVEN_IMG", run_raven),
		Tool("flye", "FLYE_IMG", run_flye),
	]
	random.shuffle(assembly_tools)  # shuffle the order

	for tool in assembly_tools:
		tool.run(sample, args, containers)
	
	#Step 3: Snippy
	snippy_jobs = [
		(sample.unicycler_assembly, sample.snippy_unicycler_dir, "ctg"),
		(sample.raven_assembly, sample.snippy_raven_dir, "ctg"),
		(sample.flye_assembly, sample.snippy_flye_dir, "ctg"),
	]

	random.shuffle(snippy_jobs)

	for input_file, outdir, mode in snippy_jobs:
		run_snippy(sample, input_file, outdir, args, containers.SNIPPY_IMG, mode)


	print(f"=== Finished sample: {sample.name} ===")



############################################
# SAMPLE TABLE PARSING
############################################

def read_sample_table(sample_table: Path, workdir: Path):
	samples = []
	with sample_table.open() as f:
		header = next(f).strip().split("\t")
		col_sample = header.index("Sample")
		col_ref = header.index("Reference_gbk")
		col_reads = header.index("Raw_reads")
		for line in f:
			row = line.strip().split("\t")
			if not row[col_sample].startswith("#"):
				samples.append(Sample(row[col_sample], Path(row[col_ref].strip()), Path(row[col_reads].strip()), workdir))
	return samples

############################################
# COUNT VARIANTS AND CREATE STATS
############################################

# noinspection PyArgumentList
def count_variants(snps_file: Path) -> pd.DataFrame:
	"""Load snps.csv and return a DataFrame with unique variants (CHROM, POS, TYPE, REF, ALT)."""
	if not snps_file.is_file():
		return pd.DataFrame(columns=KEY_COLS + ['EVIDENCE','FTYPE','STRAND','NT_POS','AA_POS','EFFECT','LOCUS_TAG','GENE','PRODUCT'])
	df = pd.read_csv(snps_file, sep=",")
	return df.drop_duplicates(subset=KEY_COLS)

def assembly_stats(assembly_file: Path):
	"""Return (number of contigs, total length) for a FASTA assembly."""
	if assembly_file.is_file() and assembly_file.stat().st_size > 0:
		n = 0
		total = 0
		for rec in SeqIO.parse(str(assembly_file), "fasta"):
			n += 1
			total += len(rec.seq)
		return n, total
	return 0, 0


def fastq_stats(fastq_file: Path):
	"""
	Return (number of reads, total nucleotides) in a FASTQ file.
	Works for .gz or plain files.
	"""
	if not fastq_file.is_file() or fastq_file.stat().st_size == 0:
		return 0, 0

	# Open gzipped or plain
	if fastq_file.suffix == ".gz":
		handle = gzip.open(fastq_file, "rt")
	else:
		handle = open(fastq_file, "r")

	n_reads, total_nt = 0, 0
	with handle:
		for rec in SeqIO.parse(handle, "fastq"):
			n_reads += 1
			total_nt += len(rec.seq)

	return n_reads, total_nt


def intersect_variants(dfs) -> pd.DataFrame:
	"""Return the intersection of multiple variant DataFrames on CHROM, POS, TYPE, REF, ALT."""
	if not dfs:
		return pd.DataFrame(columns=KEY_COLS)
	# Start with first DF keys
	base = dfs[0][KEY_COLS].drop_duplicates()
	for df in dfs[1:]:
		base = base.merge(df[KEY_COLS].drop_duplicates(), on=KEY_COLS, how='inner')
	return base.drop_duplicates()


def load_existing_stats(args):
	if args.stats.exists():
		print_info(f"Loading existing stats table: {args.stats}")
		# noinspection PyArgumentList
		return pd.read_csv(args.stats, sep="\t")
	return None

def snippy_stats(sample, existing_row = None) -> dict:
	"""Compute variant and assembly stats for a single sample and save merged variants."""
	#pdb.set_trace()
	print_info(f"Computing stats for sample: {sample.name}")
	snippy_files = {
		'unicycler': sample.snippy_unicycler_dir / "snps.csv",
		'raven': sample.snippy_raven_dir / "snps.csv",
		'flye': sample.snippy_flye_dir / "snps.csv",
		#'raw': sample.snippy_raw_dir / "snps.csv"
	}
	
	
	# ---- FASTQ + assembly stats (reuse if possible) ----
	if existing_row is not None:
		raw_reads = existing_row["raw_reads"]
		raw_nt = existing_row["raw_nt"]
		filt_reads = existing_row["filt_reads"]
		filt_nt = existing_row["filt_nt"]

		unicycler_contigs = existing_row["unicycler_contigs"]
		unicycler_total_len = existing_row["unicycler_total_len"]
		raven_contigs = existing_row["raven_contigs"]
		raven_total_len = existing_row["raven_total_len"]
		flye_contigs = existing_row["flye_contigs"]
		flye_total_len = existing_row["flye_total_len"]
	else:
		raw_reads, raw_nt = fastq_stats(sample.raw_reads)
		filt_reads, filt_nt = fastq_stats(sample.filtlong_fastq)

		unicycler_contigs, unicycler_total_len = assembly_stats(sample.unicycler_assembly)
		raven_contigs, raven_total_len = assembly_stats(sample.raven_assembly)
		flye_contigs, flye_total_len = assembly_stats(sample.flye_assembly)

	# ---- Variant handling ----
	if sample.variants.exists() and existing_row is not None:
		# Fully reuse existing values
		flye_raven_count = existing_row["flye_raven_variants"]
		raven_unicycler_count = existing_row["raven_unicycler_variants"]
		flye_unicycler_count = existing_row["flye_unicycler_variants"]
		all_three_count = existing_row["flye_raven_unicycler_variants"]

		snippy_counts = {
			"unicycler": existing_row["snippy_unicycler_variants"],
			"raven": existing_row["snippy_raven_variants"],
			"flye": existing_row["snippy_flye_variants"],
			#"raw": existing_row["snippy_raw_variants"],
		}
	else:
		# Load variants fresh
		variants = {k: count_variants(f) for k, f in snippy_files.items()}
		snippy_counts = {k: len(v) for k, v in variants.items()}
		
		# Pairwise intersections
		flye_raven = intersect_variants([variants["flye"], variants["raven"]])
		raven_unicycler = intersect_variants([variants["raven"], variants["unicycler"]])
		flye_unicycler = intersect_variants([variants["flye"], variants["unicycler"]])
		
		# Triple intersection
		all_three = intersect_variants([
			variants["flye"],
			variants["raven"],
			variants["unicycler"],
		])
		
		flye_raven_count = len(flye_raven)
		raven_unicycler_count = len(raven_unicycler)
		flye_unicycler_count = len(flye_unicycler)
		all_three_count = len(all_three)
		
		# Write merged variants ONLY if missing
		if not sample.variants.exists():
			merged = all_three.drop_duplicates()
			merged.to_csv(sample.variants, sep="\t", index=False)
	
	# Build stats dict
	stats = {
		"sample": sample.name,
		"reference": str(sample.ref_gbk).rsplit("/")[-1],

		"raw_reads": raw_reads,
		"raw_nt": raw_nt,
		"filt_reads": filt_reads,
		"filt_nt": filt_nt,

		"unicycler_contigs": unicycler_contigs,
		"unicycler_total_len": unicycler_total_len,
		"raven_contigs": raven_contigs,
		"raven_total_len": raven_total_len,
		"flye_contigs": flye_contigs,
		"flye_total_len": flye_total_len,

		"snippy_unicycler_variants": snippy_counts["unicycler"],
		"snippy_raven_variants": snippy_counts["raven"],
		"snippy_flye_variants": snippy_counts["flye"],
		#"snippy_raw_variants": snippy_counts["raw"],

		
		"flye_raven_variants": flye_raven_count,
		"raven_unicycler_variants": raven_unicycler_count,
		"flye_unicycler_variants": flye_unicycler_count,
		"flye_raven_unicycler_variants": all_three_count,
	}
	print_info(f"Finished stats for {sample.name}\n")
	return stats


# noinspection PyArgumentList
def create_stats_table(samples: list, args):
	"""Compute stats for all samples in parallel and write summary table."""
	print_info(f"Computing stats for {len(samples)} samples...")
	existing_df = None
	if args.stats.exists():
		print_info(f"Loading existing stats table: {args.stats}")
		existing_df = pd.read_csv(args.stats, sep="\t").set_index("sample")


	stats_list = []
	total = len(samples)
	completed = 0
	
	max_workers = min(args.threads, len(samples)) if not args.test else 1
	
	if args.test:
		for sample in samples:
			existing_row = (existing_df.loc[sample.name] if existing_df is not None and sample.name in existing_df.index else None)
			stats = snippy_stats(sample, existing_row)
			stats_list.append(stats)
	else:
		with ThreadPoolExecutor(max_workers=max_workers) as executor:
			future_to_sample = {}
			for sample in samples:
				existing_row = (existing_df.loc[sample.name] if existing_df is not None and sample.name in existing_df.index else None)
				future = executor.submit(snippy_stats, sample, existing_row)
				future_to_sample[future] = sample

			for future in as_completed(future_to_sample):
				sample = future_to_sample[future]
				try:
					stats = future.result()
					stats_list.append(stats)
				except Exception as e:
					print_error(f"Stats failed for sample {sample.name}: {e}")
					log_error(f"Stats failed for sample {sample.name}: {e}", args.errfile)
				finally:
					completed += 1
					print_info(f"Stats progress: {completed}/{total}")
	df = pd.DataFrame(stats_list)
	if not df.empty:
		df.sort_values("sample", inplace=True)
	else:
		print_error("No stats were generated successfully")
	
	df.to_csv(args.stats, sep="\t", index=False)
	print_info(f"Stats table written to {args.stats}")
	
	
	

def combine_all_variants(samples: list, args):
	"""Combine per-sample triple-intersection variant tables into a single file."""
	print_info("Combining all variant tables (intersection only)...")

	all_variants = []

	for sample in samples:
		if sample.variants.exists():
			df = pd.read_csv(sample.variants, sep="\t")
			if not df.empty:
				df.insert(0, "reference", sample.ref_gbk.name)
				df.insert(0, "sample", sample.name)
				all_variants.append(df)

	if args.all_variants.exists():
		print_info(f"Combined variants table already exists: {args.all_variants}")
		return

	if all_variants:
		combined = pd.concat(all_variants, ignore_index=True)
		combined.to_csv(args.all_variants, sep="\t", index=False)
		print_info(f"Combined variants table written to {args.all_variants}")
	else:
		print_info("No variants found in any sample to combine.")
	

############################################
# MAIN
############################################

def main():
	random.seed(42)
	args = parse_args()
	containers = Containers()

	# Check containers
	missing = [img for img in containers.all_images if not img.is_file()]
	if missing:
		raise FileNotFoundError("Missing container(s):\n" + "\n".join(str(m) for m in missing))
	
	if not args.dry_run:
		args.workdir.mkdir(exist_ok=True)
	samples = read_sample_table(args.sample_table, args.workdir)
	random.shuffle(samples)
	
	if args.call_variants:
		if args.test:
			for sample in samples:
				try:
					get_stats_and_variants_for_sample(sample, args, containers)
				except Exception as e:
					print_error(f"Error processing sample {sample.name}: {e}")
					log_error(f"Error processing sample {sample.name}: {e}", args.errfile)
		else:
			with ThreadPoolExecutor(max_workers=min(args.threads, len(samples))) as executor:
				futures = [executor.submit(get_stats_and_variants_for_sample, s, args, containers) for s in samples]
				for future in as_completed(futures):
					try:
						future.result()
					except Exception as e:
						print(f"Error processing sample: {e}")
	
	if args.create_stats:
		create_stats_table(samples, args)
	
	if args.create_stats or args.combine_variants:
		combine_all_variants(samples, args)

if __name__ == "__main__":
	main()
