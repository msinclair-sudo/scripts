#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 <command> [options]"
    echo
    echo "Commands:"
    echo "  exe       Extract regions from FASTA files"
    echo "  align     Align sequences"
    echo "  cluster   Cluster sequences"
    echo "  buildTree Build phylogenetic tree"
    echo
    echo "Use '$0 <command> --help' for more information on a specific command."
    exit 1
}

# Function to display usage for the 'exe' command
exe_usage() {
    echo "Usage: $0 exe [options]"
    echo
    echo "Options for 'exe' command:"
    echo "  -i      input_dir      Directory containing the input FASTA files"
    echo "  -o      output_dir     Directory to save the output files"
    echo "  -s      start_position Start position of the region to extract"
    echo "  -e      end_position   End position of the region to extract"
    echo "  --bed   bedFile.bed    BED file containing regions to extract"
    echo "  --force                Force overwrite of existing output files"
    echo "  --help                 Display this help message and exit"
    exit 1
}

# Function to display usage for the 'align' command
align_usage() {
    echo "Usage: $0 align [options]"
    echo
    echo "Options for 'align' command:"
    echo "  -i      input_dir      Directory containing the input FASTA files"
    echo "  -o      output_dir     Directory to save the output files"
    echo "  --method method        Alignment method (1, 2, 3, or 4)"
    echo "  --help                 Display this help message and exit"
    exit 1
}

# Function to display usage for the 'cluster' command
cluster_usage() {
    echo "Usage: $0 cluster [options]"
    echo
    echo "Options for 'cluster' command:"
    echo "  -i      input_dir      Directory containing the input FASTA files"
    echo "  -o      output_dir     Directory to save the output files"
    echo "  -c      threshold      Clustering threshold"
    echo "  --cd-hit-help          Display the help message for the cd-hit program"
    echo "  --help                 Display this help message and exit"
    exit 1
}

# Function to display usage for the 'buildTree' command
buildTree_usage() {
    echo "Usage: $0 buildTree [options]"
    echo
    echo "Options for 'buildTree' command:"
    echo "  -i      input_dir      Directory containing the input FASTA files"
    echo "  -o      output_dir     Directory to save the output files"
    echo "  --method method        Tree building method (1, 2, 3, or 4)"
    echo "  --help                 Display this help message and exit"
    exit 1
}

# Function to check if directory contains only .fa or .fasta files
check_fasta_files() {
    local dir="$1"
    local fa_count=$(ls "$dir"/*.fa 2>/dev/null | wc -l)
    local fasta_count=$(ls "$dir"/*.fasta 2>/dev/null | wc -l)
    local total_files=$(ls "$dir" | wc -l)

    if [ "$total_files" -ne "$fa_count" ] && [ "$total_files" -ne "$fasta_count" ]; then
        echo "Error: input_dir should contain only .fa or .fasta files, but not both."
        exit 1
    fi

    if [ "$fa_count" -gt 0 ]; then
        echo "$dir/*.fa"
    else
        echo "$dir/*.fasta"
    fi
}

# Function to execute the 'exe' command
exe_command() {
    local input_dir="$1"
    local output_dir="$2"
    local start_position="$3"
    local end_position="$4"
    local bed_file="$5"
    local force_overwrite="$6"

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Create a temporary directory within the output directory
    local temp_dir=$(mktemp -d -p "$output_dir")

    # Create the logs directory within the output directory
    local logs_dir="$output_dir/logs"
    mkdir -p "$logs_dir"

    # Create the bed file if -s and -e are provided
    if [ -z "$bed_file" ]; then
        bed_file="$output_dir/bedFile.bed"
    fi

    # Count the number of FASTA files
    local total_files=$(ls -1 "$input_dir"/*.fa "$input_dir"/*.fasta 2>/dev/null | wc -l)
    local current_file=0

    # Loop through each FASTA file in the directory
    for fasta_file in "$input_dir"/*.fa "$input_dir"/*.fasta; do
        # Increment the current file counter
        current_file=$((current_file + 1))

        # Get the base name of the file (without the directory and extension)
        local base_name=$(basename "$fasta_file" .fa)
        base_name=$(basename "$base_name" .fasta)

        # Output file name
        local output_file="$output_dir/${base_name}_trim.fa"

        # Check if the output file already exists
        if [[ -f "$output_file" ]] && [ "$force_overwrite" = false ]; then
            echo "One or more output files already exist. Use the --force argument to overwrite."
            rm -rf "$temp_dir"
            exit 1
        fi

        # Read the header line, remove whitespace, replace '|' and ':' with ',', and save to variable "header"
        local original_header=$(grep "^>" "$fasta_file" | head -n 1)
        local modified_header=$(echo "$original_header" | tr -d ' ' | tr '|' ',' | cut -d',' -f1)

        # Create a new file with the modified header in the temporary directory
        local new_fasta_file="$temp_dir/${base_name}_modified.fa"
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
            local clean_header=$(echo "$modified_header" | tr -d '>')
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
}

# Placeholder function to execute the 'align' command
align_command() {
    local input_dir="$1"
    local output_dir="$2"
    local method="$3"

    echo "Align command is not yet implemented."
    echo "Input directory: $input_dir"
    echo "Output directory: $output_dir"
    echo "Method: $method"
}

# Placeholder function to execute the 'cluster' command
cluster_command() {
    local input_dir="$1"
    local output_dir="$2"
    local threshold="$3"

    echo "Cluster command is not yet implemented."
    echo "Input directory: $input_dir"
    echo "Output directory: $output_dir"
    echo "Clustering threshold: $threshold"
}

# Placeholder function to execute the 'buildTree' command
buildTree_command() {
    local input_dir="$1"
    local output_dir="$2"
    local method="$3"

    echo "BuildTree command is not yet implemented."
    echo "Input directory: $input_dir"
    echo "Output directory: $output_dir"
    echo "Method: $method"
}

# Main script logic
if [ $# -lt 1 ]; then
    usage
fi

command=$1
shift

case $command in
    exe)
        force_overwrite=false
        while [[ "$#" -gt 0 ]]; do
            case $1 in
                -i) input_dir="$2"; shift ;;
                -o) output_dir="$2"; shift ;;
                -s) start_position="$2"; shift ;;
                -e) end_position="$2"; shift ;;
                --bed) bed_file="$2"; shift ;;
                --force) force_overwrite=true ;;
                --help) exe_usage ;;
                *) echo "Unknown parameter passed: $1"; exe_usage ;;
            esac
            shift
        done

        # Check if all required arguments are provided
        if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
            exe_usage
        fi

        # Execute the 'exe' command
        exe_command "$input_dir" "$output_dir" "$start_position" "$end_position" "$bed_file" "$force_overwrite"
        ;;
    align)
        while [[ "$#" -gt 0 ]]; do
            case $1 in
                -i) input_dir="$2"; shift ;;
                -o) output_dir="$2"; shift ;;
                --method) method="$2"; shift ;;
                --help) align_usage ;;
                *) echo "Unknown parameter passed: $1"; align_usage ;;
            esac
            shift
        done

        # Check if all required arguments are provided
        if [ -z "$input_dir" ] || [ -z "$output_dir" ] || [ -z "$method" ]; then
            align_usage
        fi

        # Execute the 'align' command
        align_command "$input_dir" "$output_dir" "$method"
        ;;
    cluster)
        while [[ "$#" -gt 0 ]]; do
            case $1 in
                -i) input_dir="$2"; shift ;;
                -o) output_dir="$2"; shift ;;
                -c) threshold="$2"; shift ;;
                --cd-hit-help) echo "CD-HIT help message"; exit 0 ;;
                --help) cluster_usage ;;
                *) echo "Unknown parameter passed: $1"; cluster_usage ;;
            esac
            shift
        done

        # Check if all required arguments are provided
        if [ -z "$input_dir" ] || [ -z "$output_dir" ] || [ -z "$threshold" ]; then
            cluster_usage
        fi

        # Execute the 'cluster' command
        cluster_command "$input_dir" "$output_dir" "$threshold"
        ;;
    buildTree)
        while [[ "$#" -gt 0 ]]; do
            case $1 in
                -i) input_dir="$2"; shift ;;
                -o) output_dir="$2"; shift ;;
                --method) method="$2"; shift ;;
                --help) buildTree_usage ;;
                *) echo "Unknown parameter passed: $1"; buildTree_usage ;;
            esac
            shift
        done

        # Check if all required arguments are provided
        if [ -z "$input_dir" ] || [ -z "$output_dir" ] || [ -z "$method" ]; then
            buildTree_usage
        fi

        # Execute the 'buildTree' command
        buildTree_command "$input_dir" "$output_dir" "$method"
        ;;
    *)
        echo "Unknown command: $command"
        usage
        ;;
esac