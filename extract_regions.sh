#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 -i input_dir -o output_dir [-s start_position -e end_position | --bed bedFile.bed]"
    echo
    echo "Options:"
    echo "  -i      input_dir      Directory containing the input FASTA files"
    echo "  -o      output_dir     Directory to save the output files"
    echo "  -s      start_position Start position of the region to extract"
    echo "  -e      end_position   End position of the region to extract"
    echo "  --bed   bedFile.bed    BED file containing regions to extract"
    echo "  --help  Help Me        Display this help message and exit"
    echo
    echo " If nothing happened, check the log file on the output dir. It might say something *shrug*"
    exit 1
}

# Check if no arguments are provided
if [ $# -eq 0 ]; then
    usage
fi

# Parse command-line arguments
while getopts ":i:o:s:e:-:" opt; do
    case $opt in
        i) input_dir="$OPTARG"
        ;;
        o) output_dir="$OPTARG"
        ;;
        s) start_position="$OPTARG"
        ;;
        e) end_position="$OPTARG"
        ;;
        -)
            case "${OPTARG}" in
                bed)
                    bed_file="${!OPTIND}"; OPTIND=$((OPTIND + 1))
                    ;;
                help)
                    usage
                    ;;
                *)
                    echo "Invalid option --${OPTARG}" >&2
                    usage
                    ;;
            esac
        ;;
        \?) echo "Invalid option -$OPTARG" >&2
            usage
        ;;
        :) echo "Option -$OPTARG requires an argument." >&2
            usage
        ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
    usage
fi

# Check for conflicting options
if [ -n "$bed_file" ] && ([ -n "$start_position" ] || [ -n "$end_position" ]); then
    echo "Error: -s and -e cannot be provided together with --bed" >&2
    usage
fi

# Check if bedtools is installed if --bed is provided
if [ -n "$bed_file" ] && ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is not installed. Please install bedtools to use the --bed option." >&2
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Create a temporary directory within the output directory
temp_dir=$(mktemp -d -p "$output_dir")

# Create the logs directory within the output directory
logs_dir="$output_dir/logs"
mkdir -p "$logs_dir"

# Create the bed file if -s and -e are provided
if [ -z "$bed_file" ]; then
    bed_file="$output_dir/bedFile.bed"
fi

# Count the number of FASTA files
total_files=$(ls -1 "$input_dir"/*.fa | wc -l)
current_file=0

# Loop through each FASTA file in the directory
for fasta_file in "$input_dir"/*.fa; do
    # Increment the current file counter
    current_file=$((current_file + 1))

    # Check if the FASTA file exists
    if [[ ! -f "$fasta_file" ]]; then
        echo "FASTA file not found: $fasta_file"
        continue
    fi

    # Get the base name of the file (without the directory and extension)
    base_name=$(basename "$fasta_file" .fa)

    # Output file name
    output_file="$output_dir/${base_name}_trim.fa"

    # Read the header line, remove whitespace, replace '|' and ':' with ',', and save to variable "header"
    original_header=$(grep "^>" "$fasta_file" | head -n 1)
    modified_header=$(echo "$original_header" | tr -d ' ' | tr '|' ',' | cut -d',' -f1)

    # Create a new file with the modified header in the temporary directory
    new_fasta_file="$temp_dir/${base_name}_modified.fa"
    echo "$modified_header" > "$new_fasta_file"
    grep -v "^>" "$fasta_file" >> "$new_fasta_file"

    # Create an index for the new FASTA file
    samtools faidx "$new_fasta_file"

    if [ -n "$bed_file" ]; then
        # Use bedtools to extract regions from the provided BED file and redirect errors to the log file
        bedtools getfasta -fi "$new_fasta_file" -bed "$bed_file" -fo "$output_file" 2>> "$logs_dir/errors.log" &
        wait $!
    else
        # Extract the specified region using the modified header (without the '>')
        clean_header=$(echo "$modified_header" | tr -d '>')
        samtools faidx "$new_fasta_file" "${clean_header}:${start_position}-${end_position}" > "$output_file"

        # Write to the bed file
        echo -e "${clean_header}\t${start_position}\t${end_position}" >> "$bed_file"
    fi

    # Display progress
    echo -ne "Processing file $current_file of $total_files\r"
done

# Remove the temporary directory
rm -rf "$temp_dir"

echo -ne "\nAll files processed.\n"