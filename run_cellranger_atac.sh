#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J cranger-atac
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --exclude me[201-248]
# #SBATCH --mail-type=NONE
# #SBATCH --mail-user=antti.kiviaho@tuni.fi

main_dir="/lustre/scratch/kiviaho/dbitseq_atac/"
cellranger_whitelist_dir="tools/cellranger-atac-1.2.0/cellranger-atac-cs/1.2.0/lib/python/barcodes/"
whitelist_name="737K-cratac-v1.txt.gz"
valid_barcodes_file="spatial_location_barcodes.csv"

sample_name="sample_A"

# Extract valid barcodes from the 'index' column and save them for cellranger to find
cut -d',' -f2 $valid_barcodes_file | tail -n +2 |  gzip > "${main_dir}${cellranger_whitelist_dir}${whitelist_name}"
#sed 's/^..//' | sed 's/..$//' |

# Run cellranger ATAC
${main_dir}tools/cellranger-atac-1.2.0/cellranger-atac count \
    --id $sample_name \
    --reference /lustre/compbio/pub/references/cellranger-atac/refdata-cellranger-atac-GRCh38-1.2.0 \
    --fastqs ${main_dir}fastqs



