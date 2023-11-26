# Code written by Danielle Rivera
# Purpose: extracting repeat elements fr a gff file generated using RepeatMasker

# Code was created in part with assistance from OpenAI:
# Author: OpenAI
# Title: Custom Python Script Assistance
# Year: 2023
# URL: https://www.openai.com/


import csv

def convert_gff_to_csv(gff_file, csv_file):
    with open(gff_file, 'r') as gff, open(csv_file, 'w', newline='') as csv_output:
        csv_writer = csv.writer(csv_output, delimiter='\t')
        csv_writer.writerow(['Chromosome', 'Start', 'End'])  # Header

        for line in gff:
            # Skip comments and header lines
            if line.startswith('#') or line.strip() == '':
                continue

            fields = line.strip().split('\t')
            chromosome = fields[0]
            start = fields[3]
            end = fields[4]

            csv_writer.writerow([chromosome, start, end])

# Specify your RepeatMasker GFF file and the desired CSV output file
gff_file = "Snitidus_v0.9.fasta.repeats_out.gff"
csv_output_file = "Snitidus_repeatmasker_output.csv"

# Convert GFF to CSV
convert_gff_to_csv(gff_file, csv_output_file)

print(f"Conversion complete. CSV file saved to {csv_output_file}")
