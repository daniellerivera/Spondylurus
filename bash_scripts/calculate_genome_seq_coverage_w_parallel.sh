#!/bin/bash

# Create the output CSV file
echo "Sample,Total Bases,Genome Coverage" > sample_sequencing_coverage.csv

# Function to calculate total bases and genome coverage for a folder
calculate_coverage() {
    folder=$1
    folder_name=$(basename "$folder")

    # Check if there are fq.gz files in the folder
    if ls "$folder"/*.fq.gz 1> /dev/null 2>&1; then
        # Count the total number of bases in all fq.gz files
        total_bases=$(zcat "$folder"/*.fq.gz | awk 'NR%4 == 2 {sum += length($1)} END {print sum}')

        # Calculate genome coverage
        genome_length=1390000000
        genome_coverage=$(bc <<< "scale=10; $total_bases / $genome_length")

        # Print results to CSV
        echo "$folder_name,$total_bases,$genome_coverage" >> sample_sequencing_coverage.csv
    else
        echo "No fq.gz files found in $folder_name"
    fi
}

# Export the function for parallel execution
export -f calculate_coverage

# Find all subdirectories (excluding the current directory) and run the function in parallel
find . -maxdepth 1 -mindepth 1 -type d -exec bash -c 'calculate_coverage "$0"' {} \;
