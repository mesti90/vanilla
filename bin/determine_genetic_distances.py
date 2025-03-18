import pandas as pd
from scipy.spatial.distance import pdist, squareform
import argparse
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch

def load_variants(file_path):
	"""Load the TSV file containing variant data."""
	return pd.read_csv(file_path, sep="\t")

def load_all_samples(sample_list_file):
	"""Load sample names and reference genome information from a file."""
	sample_df = pd.read_csv(sample_list_file, sep="\t")
	return sample_df["Sample"].tolist(), sample_df["Reference_name"].unique()[0]  # Extract reference genome

def filter_variants(df):
	"""Filter SNPs and Indels from the dataset."""
	del_ins_types = ["snp", "del", "ins"]
	return df[df["TYPE"].isin(del_ins_types)]

def create_variant_matrix(variant_df, all_samples, reference_name):
	"""Create a pivot table with Samples as rows and variant positions as columns."""
	variant_matrix = variant_df.pivot_table(index="Sample", columns=["CHROM", "POS"], values="ALT", aggfunc="first", fill_value=".")
	
	# Ensure all samples are included in the matrix, even if they had no variants
	for sample in all_samples + [reference_name]:  # Include reference
		if sample not in variant_matrix.index:
			variant_matrix.loc[sample] = "."  # No variants for this sample
	
	return variant_matrix

def compute_distance_matrix(variant_matrix):
	"""Compute pairwise Hamming distances from the variant matrix."""
	binary_matrix = (variant_matrix != ".").astype(int)
	
	# Compute distances
	if binary_matrix.shape[0] > 1:  # Ensure we have more than one sample
		distance_matrix = pd.DataFrame(
			squareform(pdist(binary_matrix, metric="hamming") * binary_matrix.shape[1]),
			index=binary_matrix.index,
			columns=binary_matrix.index
		)
	else:
		# If only one sample, return a zero-distance matrix
		distance_matrix = pd.DataFrame(0, index=binary_matrix.index, columns=binary_matrix.index)
	
	return distance_matrix

def build_tree(distance_matrix, tree_output):
	"""Generate a hierarchical clustering tree from the distance matrix."""
	if distance_matrix.shape[0] > 1:  # Only build a tree if there are multiple samples
		condensed_distances = pdist(distance_matrix, metric="euclidean")
		linkage_matrix = sch.linkage(condensed_distances, method='average')
		plt.figure(figsize=(10, 6))
		sch.dendrogram(linkage_matrix, labels=distance_matrix.index, leaf_rotation=60)
		plt.title("Hierarchical Clustering Dendrogram (SNP Distance)")
		plt.xlabel("Samples")
		plt.ylabel("Genetic Distance")
		plt.savefig(tree_output)
		plt.show()
	else:
		print("Not enough samples for tree construction.")

def main():
	parser = argparse.ArgumentParser(description="Compute pairwise SNP and Indel distances from a TSV file and generate a phylogenetic tree.")
	parser.add_argument("-i", "--input_file", help="Path to the input TSV file", default="Variant_results/pseudomonas.variants.20241125.tsv")
	parser.add_argument("-s", "--sample_list", help="Path to the file containing a list of all sample names", default="config/vanilla.samples.20241125.tsv")
	parser.add_argument("-o", "--output_file", help="Path to save the output distance matrix", default="Variant_results/pseudomonas.distances.20250318.tsv")
	parser.add_argument("-t", "--tree_output", help="Path to save the output dendrogram image", default="Variant_results/pseudomonas.tree.20250318.pdf")
	args = parser.parse_args()

	df = load_variants(args.input_file)
	all_samples, reference_name = load_all_samples(args.sample_list)
	variant_df = filter_variants(df)
	variant_matrix = create_variant_matrix(variant_df, all_samples, reference_name)
	distance_matrix = compute_distance_matrix(variant_matrix)
	
	# Save the distance matrix to a file
	distance_matrix.to_csv(args.output_file, sep="\t")
	
	# Build and save the tree
	build_tree(distance_matrix, args.tree_output)
	
	# Display the matrix
	print(distance_matrix)

if __name__ == "__main__":
	main()
