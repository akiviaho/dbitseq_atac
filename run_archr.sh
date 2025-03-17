#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J archR
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
# #SBATCH --exclude me[201-248]
# #SBATCH --mail-type=NONE
# #SBATCH --mail-user=antti.kiviaho@tuni.fi

main_dir="/lustre/scratch/kiviaho/dbitseq_atac/"
sample_name="sample_A"

module load anaconda
source activate archr

Rscript archr_preprocessing.R

