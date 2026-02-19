# Code written by Danielle Rivera
# Purpose: creating a file that measures GC-content across a genome, to be used to create a Circos Plot

# Code was created in part with assistance from OpenAI:
# Author: OpenAI
# Title: Custom Python Script Assistance
# Year: 2023
# URL: https://www.openai.com/

from Bio import SeqIO

def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    gc_content = (gc_count / total_bases) * 100
    return gc_content

def bin_gc_content(chromosome_name, reference_genome, bin_size):
    gc_content_values = []
    genome_length = len(reference_genome)

    for start in range(0, genome_length, bin_size):
        end = min(start + bin_size, genome_length)
        window_sequence = reference_genome[start:end]
        gc_content = calculate_gc_content(window_sequence)
        gc_content_values.append((chromosome_name, start, end, gc_content))

    return gc_content_values

# Specify your genome file name
genome_file = "Snitidus_ref_1.0.fasta"

# Set the bin size (500 kilobases)
bin_size = 500 * 10**3

# Read the reference genome
gc_content_data = []
with open(genome_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        chromosome_name = record.id
        reference_genome = str(record.seq)
        gc_content_values = bin_gc_content(chromosome_name, reference_genome, bin_size)
        gc_content_data.extend(gc_content_values)

# Append GC content values to an output file
output_file = "gc_content_binned_with_chromosome.txt"
with open(output_file, "w") as f:
    f.write("Chromosome\tStart\tEnd\tGC_Content\n")
    for chromosome_name, start, end, gc_content in gc_content_data:
        f.write(f"{chromosome_name}\t{start}\t{end}\t{gc_content}\n")

print(f"GC content values have been written to {output_file}")