#!/bin/bash

module load compbio/samtools

threads=8
sample_name="sample_A"

sort --parallel $threads -k 1,1 -k2,2n -o ${sample_name}_Sorted.bed ${sample_name}_alignment.bed;
bgzip -@ $threads ${sample_name}_Sorted.bed;
tabix -p bed ${sample_name}_Sorted.bed.gz
