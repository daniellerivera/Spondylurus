# Code written by Danielle Rivera
# Purpose: splitting a large csv file into multiple files for downstream analyses

# Code was created in part with assistance from OpenAI:
# Author: OpenAI
# Title: Custom Python Script Assistance
# Year: 2023
# URL: https://www.openai.com/

import csv
import os

def split_csv(input_csv, output_folder, max_lines=2500000):
    with open(input_csv, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        header = next(csv_reader)

        current_chunk = 1
        current_lines = 0

        while True:
            output_csv = os.path.join(output_folder, f"Snitidus_repeatmasker_out_{current_chunk}.csv")

            with open(output_csv, 'w', newline='') as chunk_file:
                csv_writer = csv.writer(chunk_file, delimiter=',')
                csv_writer.writerow(header)

                for _ in range(max_lines):
                    try:
                        row = next(csv_reader)
                    except StopIteration:
                        break
                    csv_writer.writerow(row)
                    current_lines += 1

            if current_lines == 0:
                # No more rows to process
                break

            current_chunk += 1
            current_lines = 0

# Specify your input CSV file and output folder
input_csv = "Snitidus_repeatmasker_out.csv"
output_folder = "output_chunks"

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Split the CSV file into chunks
split_csv(input_csv, output_folder)

print(f"CSV file split into chunks in the {output_folder} directory.")
