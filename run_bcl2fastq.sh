#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J bcl2fastq
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
# #SBATCH --mail-type=NONE
# #SBATCH --mail-user=antti.kiviaho@tuni.fi

# Default values for variables
bclDir=""
outdir=""
format_names_for_cellranger=true

# Function to display usage
usage() {
    echo "Usage: $0 --bclDir <bclDir> --outdir <outdir> --sample_name <sample_name> --format_names_for_cellranger <bool>"
    exit 1
}

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bclDir) bclDir="$2"; shift ;;
        --outdir) outdir="$2"; shift ;;
        --sample_name) sample_name="$2"; shift ;;
        --format_names_for_cellranger) format_names_for_cellranger=$2; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if both parameters are provided
if [ -z "$bclDir" ] || [ -z "$outdir" ] || [ -z "$sample_name" ]; then
    usage
fi

/lustre/scratch/kiviaho/dbitseq_atac/tools/bin/bcl2fastq \
-R "$bclDir" \
-o "$outdir" \
--create-fastq-for-index-reads \
--no-lane-splitting \
-r 32 -p 32 -w 32 \
--mask-short-adapter-reads 0

# Rename the resulting files
for filename in "$outdir"/*; do
    # Check if the filename contains 'Undetermined_S0'
    if [[ "$filename" == *"Undetermined_S0"* ]]; then
        # Create the new filename by replacing 'Undetermined_S0' with the value of sample_name
        new_filename="${filename//Undetermined_S0/$sample_name}"
        
        # Additional replacements if format_names_for_cellranger is true
        if [[ "$format_names_for_cellranger" == true ]]; then
            new_filename="${new_filename//R1_001.fastq.gz/S1_L001_R1_001.fastq.gz}"
            new_filename="${new_filename//R2_001.fastq.gz/S1_L001_R3_001.fastq.gz}"
            new_filename="${new_filename//I1_001.fastq.gz/untrimmed_S1_L001_R2_001.fastq.gz}"
        fi
        # Rename the file
        mv "$filename" "$new_filename"
        echo "Renamed: $filename -> $new_filename"
    fi
done

