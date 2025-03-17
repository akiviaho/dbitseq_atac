#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J chromap_bam
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
# #SBATCH --exclude me[201-248]
# #SBATCH --mail-type=NONE
# #SBATCH --mail-user=antti.kiviaho@tuni.fi

main_dir="/lustre/scratch/kiviaho/dbitseq_atac/"
valid_barcodes_file="spatial_location_barcodes.csv"
sample_name="sample_A"

module load anaconda
module load compbio/samtools
source activate chromap

# Run chromap
# Note that the files have been named according to the cellranger 
# convention: read1 = R1, read2 = R3, barcodes = R2
chromap \
    --preset atac \
    --bc-error-threshold 3 \
    --ref /lustre/compbio/pub/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
    --index ${main_dir}chromap_index \
    --read1 ${main_dir}fastqs/${sample_name}_S1_L001_R1_001.fastq.gz \
    --read2 ${main_dir}fastqs/${sample_name}_S1_L001_R3_001.fastq.gz \
    --barcode ${main_dir}fastqs/${sample_name}_S1_L001_R2_001.fastq.gz \
    --barcode-whitelist ${main_dir}${sample_name}_barcode_whitelist.txt \
    --SAM -o /dev/stdout | samtools view -bS - | samtools sort -o ${sample_name}_sorted.bam

samtools index ${sample_name}_sorted.bam