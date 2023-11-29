import csv
import pandas as pd
from Bio import SeqIO
import fileinput
import os
from Bio.SeqUtils import GC

def import_data(file_path):
    return pd.read_csv(file_path)
    
# Specify your input files
genome_file = "SponNit_1.0.fasta"
gff_file_path = "Snitidus.gff"
scaflens_file_path = "SponNit_1.0.scaflens"
repeatmasker_gff_path = "Snitidus_v0.9.fasta.repeats_out.gff"

# Define all Output filenames
gc_output_file = "gc_content_binned_with_chromosome.txt"
variants_output_path = "variants_extracted_from_braker.csv"
chromosomes_output_path = "chromosomes.csv"
repeats_output_path = "repeats_extracted_from_repeatmasker.csv"

# Modify files with Chromosome names
chr_rename = "_Sni"

# List of files to modify
files_to_modify = [genome_file, gff_file_path, scaflens_file_path, repeatmasker_gff_path]

# Perform the renaming in each file
for file_path in files_to_modify:
    # Check if the file exists before modifying
    if os.path.exists(file_path):
        with fileinput.FileInput(file_path, inplace=True, backup=".bak") as file:
            for line in file:
                print(line.replace(chr_rename, ""), end="")
        print(f"Chromosome names modified in {file_path}")
    else:
        print(f"File {file_path} not found.")

print("Chromosome names modified successfully in all specified files.")

# Read scaflens file
scaflens_data = pd.read_csv(scaflens_file_path, sep='\t', header=None)
scaflens_data.columns = ["Chromosome", "End", "Size", "Percentage", "Order", "Details", "Telomeres"]

# Filter and format chromosomes
chromosomes_data = scaflens_data[scaflens_data["Chromosome"].str.startswith("Chr")].copy()
chromosomes_data["Start"] = 1

# Write chromosomes to a CSV file
chromosomes_data[["Chromosome", "Start", "End"]].to_csv(chromosomes_output_path, index=False)




            
def create_bins(chromosome_sizes, bin_size_mb=20):
    bin_edges = {}

    for chromosome, (start, end) in chromosome_sizes.items():
        num_bins = int((end - start) / (bin_size_mb * 1000000)) + 1
        bin_edges[chromosome] = [(i * bin_size_mb * 1000000 + start, (i + 1) * bin_size_mb * 1000000 + start) for i in range(num_bins)]

    return bin_edges
    
# Read chromosomes data
chromosomes_data = pd.read_csv(chromosomes_output_path)

# Create bins for each chromosome using information from chromosomes_data
chromosome_sizes = {row["Chromosome"]: (1, row["End"]) for _, row in chromosomes_data.iterrows()}
bin_size_mb = 1
bin_edges = create_bins(chromosome_sizes, bin_size_mb)

def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    gc_content = (gc_count / total_bases) * 100
    return gc_content

def bin_gc_content(chromosome_name, reference_genome, gc_bin_size):
    gc_content_values = []
    genome_length = len(reference_genome)

    for start in range(1, genome_length, gc_bin_size):
        end = min(start + gc_bin_size, genome_length)
        window_sequence = reference_genome[start:end]
        gc_content = calculate_gc_content(window_sequence)
        gc_content_values.append((chromosome_name, start, end, gc_content))

    return gc_content_values

# Set the bin size (500 kilobases)
gc_bin_size = 500 * 10**3

# Read the reference genome
gc_content_data = []
with open(genome_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        chromosome_name = record.id
        reference_genome = str(record.seq)
        gc_content_values = bin_gc_content(chromosome_name, reference_genome, gc_bin_size)
        gc_content_data.extend(gc_content_values)

# Append GC content values to an output file
with open(gc_output_file, "w") as f:
    f.write("Chromosome\tStart\tEnd\tGC_Content\n")
    for chromosome_name, start, end, gc_content in gc_content_data:
        f.write(f"{chromosome_name}\t{start}\t{end}\t{gc_content}\n")

print(f"GC content values have been written to {gc_output_file}")

# Read GFF3 file with variants/feature
gff_data = pd.read_csv(gff_file_path, sep='\t', header=None)
gff_data.columns = ["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"]

# Extract relevant columns
variants_data = gff_data[["Chromosome", "Start", "End", "Feature"]]

# Write variants to a CSV file
variants_data.to_csv(variants_output_path, index=False)

# Create separate files for each variant
variant_groups = variants_data.groupby("Feature")
for feature, group_data in variant_groups:
    variant_output_path = f"{feature}_extracted_from_braker.csv"
    group_data.to_csv(variant_output_path, index=False)

# Read the GFF file from Repeatmasker
repeatmasker_data = pd.read_csv(repeatmasker_gff_path, sep='\t', header=None, comment='#', names=["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"])

# Extract features from Repeatmasker
repeats_data = repeatmasker_data[["Chromosome", "Start", "End"]].copy()
repeats_data["Feature"] = "repeat"

# Write the extracted Repeatmasker features to a CSV file
repeats_data.to_csv(repeats_output_path, index=False)

# Print success messages
print(f"Variants data saved to {variants_output_path}")
print(f"Chromosomes data saved to {chromosomes_output_path}")
print(f"Repeat data saved to {repeats_output_path}")

# Combine variants and repeats data
all_variants_data = pd.concat([variants_data, repeats_data])

# Save combined data to "all_variants_in_genome.csv"
all_variants_output_path = "all_variants_in_genome.csv"
all_variants_data.to_csv(all_variants_output_path, index=False)

print(f"All variants data saved to {all_variants_output_path}")

# Specify your input files
all_variants_output_path = "all_variants_in_genome.csv"

all_variants_data = pd.read_csv(all_variants_output_path)

# Remove rows that begin with "scaf" in the "Chromosome" column
all_variants_data = all_variants_data[~all_variants_data["Chromosome"].str.startswith("scaf")]

# Save the modified data to a new CSV file
filtered_variants_output_path = "filtered_variants_in_genome.csv"
all_variants_data.to_csv(filtered_variants_output_path, index=False)

print(f"Filtered variants data saved to {filtered_variants_output_path}")

# Specify your input files
def calculate_genome_length(genome_file):
    genome_length = 0
    with open(genome_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            genome_length += len(record.seq)
    return genome_length

def create_bins(chromosome_sizes, bin_size_mb=10):
    bin_edges = {}
    for chromosome, (start, end) in chromosome_sizes.items():
        num_bins = int((end - start) / (bin_size_mb * 1000000)) + 1
        bin_edges[chromosome] = [(i * bin_size_mb * 1000000 + start, (i + 1) * bin_size_mb * 1000000 + start) for i in range(num_bins)]
    return bin_edges

def count_variants_in_bins(data, bin_edges):
    bin_counts = {feature: [0] * len(edges) for feature in data["Feature"].unique() for edges in bin_edges.values()}

    for chromosome, edges in bin_edges.items():
        for i, (bin_start, bin_end) in enumerate(edges):
            bin_data = data[(data["Chromosome"] == chromosome) & (data["Start"] >= bin_start) & (data["Start"] <= bin_end)]
            for feature in bin_counts.keys():
                bin_counts[feature][i] += (bin_data["Feature"] == feature).sum()

    return bin_counts

# Read chromosomes data
chromosomes_data = pd.read_csv(chromosomes_output_path)
all_results_data = pd.DataFrame()

# Create bins for each chromosome using information from chromosomes_data
chromosome_sizes = {row["Chromosome"]: (1, row["End"]) for _, row in chromosomes_data.iterrows()}
bin_size_mb = 10  # Adjust as needed
bin_edges = create_bins(chromosome_sizes, bin_size_mb)

# Create an empty DataFrame to store min and max values for each variant
variant_min_max_data = pd.DataFrame(columns=["Variant", "Min", "Max"])

# Count variants in bins for each chromosome
for chromosome, edges in bin_edges.items():
    # Read the filtered variants data for the current chromosome
    filtered_variants_output_path = f"filtered_variants_in_genome.csv"
    filtered_variants_data = pd.read_csv(filtered_variants_output_path)

    # Count variants in bins
    bin_counts = count_variants_in_bins(filtered_variants_data, {chromosome: edges})

    # Create a DataFrame for the result
    result_data = pd.DataFrame(bin_counts)

    # Add "Chromosome", "Start", "End" columns to the result
    result_data["Chromosome"] = [chromosome] * len(result_data)
    result_data["Start"] = [start for start, end in edges]
    result_data["End"] = [end for start, end in edges]

    # Append the result for the current chromosome to the accumulated DataFrame
    all_results_data = pd.concat([all_results_data, result_data], ignore_index=True)
	
	# Update the min and max values for each variant
    for feature in bin_counts.keys():
    	if feature not in variant_min_max_data["Variant"].values:
    	    # If the variant is not in the DataFrame, add a new row
    	    variant_min_max_data = variant_min_max_data.append({"Variant": feature, "Min": result_data[feature].min(), "Max": result_data[feature].max()}, ignore_index=True)
    	else:
    	    # If the variant is already in the DataFrame, update min and max values
    	    variant_min_max_data.loc[variant_min_max_data["Variant"] == feature, ["Min", "Max"]] = [
    	        min(variant_min_max_data.loc[variant_min_max_data["Variant"] == feature, "Min"].min(), result_data[feature].min()),
    	        max(variant_min_max_data.loc[variant_min_max_data["Variant"] == feature, "Max"].max(), result_data[feature].max())
    	    ]
        
# Reorder columns
all_results_data = all_results_data[["Chromosome", "Start", "End"] + list(bin_counts.keys())]

# Save to a single CSV file
all_results_data.to_csv("variant_counts_in_bins.csv", index=False)
print("Results for all chromosomes saved to variant_counts_in_bins.csv")

# Save individual raw data count CSV files for each variant
for feature in bin_counts.keys():
    raw_data_path = f"raw_data_count_{feature}_in_genome_all_variants.csv"
    feature_raw_data = all_results_data[["Chromosome", "Start", "End", feature]]
    feature_raw_data.to_csv(raw_data_path, index=False)
    print(f"Raw data count for {feature} in all variants saved to {raw_data_path}")
    
# Save variant_min_max_data to "variant_min_max.csv"
variant_min_max_data.to_csv("variant_min_max.csv", index=False)
print("Min and Max values for each variant saved to variant_min_max.csv")


# Specify the output directory for normalized variant count files
normalized_output_directory = "normalized_variant_counts"

# Create the output directory if it doesn't exist
os.makedirs(normalized_output_directory, exist_ok=True)

# List all raw_data_path files
raw_data_files = [
    "raw_data_count_gene_in_genome_all_variants.csv",
    "raw_data_count_mRNA_in_genome_all_variants.csv",
    "raw_data_count_exon_in_genome_all_variants.csv",
    "raw_data_count_CDS_in_genome_all_variants.csv",
    "raw_data_count_repeat_in_genome_all_variants.csv"
    # Add more files as needed
]

# Initialize a DataFrame to store normalized variant counts
merged_normalized_data = None


# Normalize the last column of each file and save each normalized file
for raw_data_path in raw_data_files:
    # Load the raw data
    raw_data = pd.read_csv(raw_data_path)

    # Extract the last column name
    last_column = raw_data.columns[-1]

    # Extract the count column for the last variant
    count_column = raw_data[last_column]

    # Convert count_column to numeric (ignore errors for non-numeric values)
    count_column = pd.to_numeric(count_column, errors='coerce')

    # Normalize the count column between 0 and 1
    min_value = count_column.min()
    max_value = count_column.max()

    if min_value != max_value:  # Avoid division by zero
        normalized_data = raw_data.copy()
        normalized_data[last_column] = ((count_column - min_value) / (max_value - min_value)).round(3)
    else:
        normalized_data = raw_data.copy()
        normalized_data[last_column] = 0

    # Save normalized count column to a new file
    normalized_file_path = os.path.join(normalized_output_directory, f"{os.path.basename(raw_data_path)}_normalized.csv")
    normalized_data.to_csv(normalized_file_path, index=False)
    print(f"Normalized count column in {os.path.basename(raw_data_path)} saved to {normalized_file_path}")

    # Merge the normalized count column to the overall DataFrame
    if merged_normalized_data is None:
        merged_normalized_data = normalized_data[["Chromosome", "Start", "End", last_column]]
    else:
        merged_normalized_data = pd.merge(merged_normalized_data, normalized_data[["Chromosome", "Start", "End", last_column]],
                                          on=["Chromosome", "Start", "End"], how="outer")

# Save merged normalized variant counts to a new CSV file
merged_normalized_output_path = "merged_normalized_variant_counts.csv"
merged_normalized_data.to_csv(merged_normalized_output_path, index=False)

print(f"Merged normalized variant counts saved to {merged_normalized_output_path}")

