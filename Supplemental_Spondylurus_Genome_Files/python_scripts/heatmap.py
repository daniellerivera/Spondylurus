# Code written by Danielle Rivera
# Purpose: creating a heatmap file for genome variants, to be used to create a Circos Plot

# Code was created in part with assistance from OpenAI:
# Author: OpenAI
# Title: Custom Python Script Assistance
# Year: 2023
# URL: https://www.openai.com/


import csv
import pandas as pd

def import_data(file_path):
    return pd.read_csv(file_path)

def min_max_normalize(counts):
    min_val = min(counts)
    max_val = max(counts)

    normalized_counts = [(2 * (x - min_val) / (max_val - min_val)) - 1 for x in counts]

    return normalized_counts

def create_bins(chromosome_sizes, bin_size_mb=20):
    bin_edges = {}

    for chromosome, (start, end) in chromosome_sizes.items():
        num_bins = int((end - start) / (bin_size_mb * 1000000)) + 1
        bin_edges[chromosome] = [(i * bin_size_mb * 1000000 + start, (i + 1) * bin_size_mb * 1000000 + start) for i in range(num_bins)]

    return bin_edges

def count_variants_in_bins(data, bin_edges):
    bin_counts = {variant: {chromosome: [0] * len(edges) for chromosome, edges in bin_edges.items()} for variant in ["Exon", "mRNA", "Repeat"]}

    for chromosome, edges in bin_edges.items():
        chromosome_data = data[data["Chromosome"] == chromosome]

        for i, (bin_start, bin_end) in enumerate(edges):
            bin_data = chromosome_data[(chromosome_data["Start"] >= bin_start) & (chromosome_data["Start"] <= bin_end)]

            for variant in ["Exon", "mRNA", "Repeat"]:
                bin_counts[variant][chromosome][i] = (bin_data["Variant"] == variant).sum()

    return bin_counts

def normalize_counts(bin_counts):
    for variant, variant_counts in bin_counts.items():
        for chromosome, counts in variant_counts.items():
            bin_counts[variant][chromosome] = min_max_normalize(counts)

    return bin_counts

def create_heatmap(bin_edges, bin_counts, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter='\t')
        csv_writer.writerow(["Chromosome", "Bin", "Start", "End", "Exon", "mRNA", "Repeat"])

        for chromosome, edges in bin_edges.items():
            for i, (bin_start, bin_end) in enumerate(edges):
                row = [
                    chromosome,
                    i + 1,
                    bin_start,
                    bin_end,
                    bin_counts["Exon"][chromosome][i],
                    bin_counts["mRNA"][chromosome][i],
                    bin_counts["Repeat"][chromosome][i]
                ]
                csv_writer.writerow(row)

# Specify your input files
gff_file_path = "Snitidus.gff.csv"
repeatmasker_file1_path = "output_chunks/Snitidus_repeatmasker_out_1.csv"
repeatmasker_file2_path = "output_chunks/Snitidus_repeatmasker_out_2.csv"

# Import data from CSV files
gff_data = import_data(gff_file_path)
repeatmasker_data1 = import_data(repeatmasker_file1_path)
repeatmasker_data2 = import_data(repeatmasker_file2_path)

# Combine data from different files
combined_data = pd.concat([gff_data, repeatmasker_data1, repeatmasker_data2])

# Define chromosome sizes
chromosome_sizes = {
    "Chr1": (1, 299172356),
    "Chr2": (1, 255164795),
    "Chr3": (1, 208671581),
    "Chr4": (1, 199700853),
    "Chr5": (1, 101272343),
    "Chr6": (1, 68056104),
    "Chr7": (1, 63778719),
    "Chr8": (1, 53019942),
    "Chr9": (1, 37978314),
    "Chr10": (1, 23371246),
    "Chr11": (1, 21874651),
    "Chr12": (1, 17641550),
    "Chr13": (1, 17051481),
    "Chr14": (1, 10970200),
    "Chr15": (1, 10901850),
    "Chr16": (1, 8873110),
}

# Create bins for each chromosome
bin_size_mb = 1
bin_edges = create_bins(chromosome_sizes, bin_size_mb)

# Count variants in each bin
bin_counts = count_variants_in_bins(combined_data, bin_edges)

# Create a heatmap file of the raw data
output_raw_heatmap_file = "heatmap_output_rawdata.txt"
create_heatmap(bin_edges, bin_counts, output_raw_heatmap_file)

# Append min and max of each variant to the raw data heatmap file
with open(output_raw_heatmap_file, 'a', newline='') as csvfile:
    csv_writer = csv.writer(csvfile, delimiter='\t')

    for variant in ["Exon", "mRNA", "Repeat"]:
        min_value = min(min(bin_counts[variant][chromosome]) for chromosome in bin_counts[variant].keys())
        max_value = max(max(bin_counts[variant][chromosome]) for chromosome in bin_counts[variant].keys())

        csv_writer.writerow([f"Min_{variant}", min_value])
        csv_writer.writerow([f"Max_{variant}", max_value])

print(f"Heatmap raw data saved to {output_raw_heatmap_file}")

# Normalize variant counts
for variant in ["Exon", "mRNA", "Repeat"]:
    for chromosome in bin_counts[variant].keys():
        counts = [bin_counts[variant][chromosome][i] for i in range(len(bin_edges[chromosome]))]
        bin_counts[variant][chromosome] = min_max_normalize(counts)

# Create a heatmap file
output_heatmap_file = "heatmap_output.txt"
create_heatmap(bin_edges, bin_counts, output_heatmap_file)

print(f"Heatmap data saved to {output_heatmap_file}")
