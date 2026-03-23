#!/usr/bin/env python3

import pandas as pd
import argparse


def load_table(path):
	"""Load variant table."""
	print("Loading variant table...")
	df = pd.read_csv(path, sep="\t")
	print(f"Loaded {len(df)} rows")
	return df


def build_variant_id(df):
	"""Create a unique identifier for variants."""
	df["VARIANT_ID"] = df[["CHROM","POS","TYPE","REF","ALT"]].astype(str).agg(":".join, axis=1)
	return df


def get_reference_variants(df, reference):
	"""Return reference dataframe and variant set."""
	ref_df = df[df["Sample"] == reference]

	if ref_df.empty:
		raise ValueError(f"Reference sample '{reference}' not found")

	ref_variants = set(ref_df["VARIANT_ID"])

	print(f"Reference sample '{reference}' has {len(ref_variants)} variants")

	return ref_df, ref_variants


def find_differences(df, ref_df, ref_variants, reference, samples=None):
	"""Find variants differing from the reference."""

	if samples:
		df = df[df["Sample"].isin(samples + [reference])]

	print("Comparing samples...")

	results = []

	sample_list = df["Sample"].unique()

	for sample in sample_list:

		if sample == reference:
			continue

		sample_df = df[df["Sample"] == sample]
		sample_variants = set(sample_df["VARIANT_ID"])

		novel = sample_variants - ref_variants
		missing = ref_variants - sample_variants

		if novel:
			tmp = sample_df[sample_df["VARIANT_ID"].isin(novel)].copy()
			tmp["STATUS"] = "novel"
			results.append(tmp)

		if missing:
			tmp = ref_df[ref_df["VARIANT_ID"].isin(missing)].copy()
			tmp["Sample"] = sample
			tmp["STATUS"] = "missing"
			results.append(tmp)

	if results:
		diff_df = pd.concat(results)
	else:
		diff_df = pd.DataFrame()

	print(f"Found {len(diff_df)} differing variants")

	return diff_df


def write_output(df, output):
	"""Write results to TSV."""
	print(f"Writing output to {output}")

	if "VARIANT_ID" in df.columns:
		df = df.drop(columns=["VARIANT_ID"])

	df.to_csv(output, sep="\t", index=False)


def parse_args():
	parser = argparse.ArgumentParser(description="Compare variants across samples using pandas")
	parser.add_argument("-t", "--table", help="Variant table (TSV)")
	parser.add_argument("-r","--reference", required=True, help="Reference sample")
	parser.add_argument("-s","--samples", nargs="*", help="Samples to compare")
	parser.add_argument("-o", "--output", default="variant_differences.tsv")
	return parser.parse_args()


def main():

	args = parse_args()

	df = load_table(args.table)

	df = build_variant_id(df)

	ref_df, ref_variants = get_reference_variants(df, args.reference)

	diff_df = find_differences(df, ref_df, ref_variants, args.reference, args.samples)

	write_output(diff_df, args.output)

	print("Done.")


if __name__ == "__main__":
	main()