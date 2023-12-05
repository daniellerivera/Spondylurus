#!/bin/bash

# Create the output CSV file
echo "Sample,Total Bases,Genome Coverage" > sample_sequencing_coverage.csv

# Loop through each subdirectory in the current directory
for folder in */; do
    # Navigate into the folder
    cd "$folder" || exit

    # Get the folder name
    folder_name=$(basename "$folder")

    # Count the total number of bases in all fq.gz files
    total_bases=$(zcat *.fq.gz | awk 'NR%4 == 2 {sum += length($1)} END {print sum}')

    # Calculate genome coverage
    genome_length=1390000000
    genome_coverage=$(bc <<< "scale=10; $total_bases / $genome_length")

    # Print results to CSV
    echo "$folder_name,$total_bases,$genome_coverage" >> ../sample_sequencing_coverage.csv

    # Move back to the parent directory
    cd ..
done
