#!/bin/bash
# This script will calculate genome coverage sequenced (e.g. 10X genome coverage) by counting the number
# of bases sequenced (in .fq.gz files) and dividing that by the expected genome length. 
# NOTE 1: this assumes your folder structure as follows:
# (1) script in directory
# (2) subfolders with fq.gz files in directory with preferred naming scheme
# NOTE 2: this will not read fq.gz files within the parent directory.
# NOTE 3: remember to change the genome length to your particular genome
# NOTE 4: change your file extension as necessary (e.g. *.fq.gz, *.fastq.gz); I haven't tried it but I 
# assume you could also alter line 26 to "cat *.fq" or "cat *.fastq" for uncompressed files


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
