# Date: 5.9.2024
# Author: Antti Kiviaho
# antti.kiviaho@tuni.fi

'''
This script is for modifying a table of base sequences to create a N by M grid as per the principles of DBiT-seq.
The barcodes should be provided in an excel table with columns 'Barcode 1' and 'Barcode 2'. The number of bases
clipped from the start and the end of the sequences is defined with the parameters n_clip_start and n_clip_end
'''

import pandas as pd
import os
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Create a DBiT-seq location grid out of raw barcodes provided in an excel sheet.\nThe resulting combinations are saved into a csv.')
parser.add_argument('--xlsx_sheet', type=str, default=None, help='Path to the Excel sheet. Should include columns "Barcode 1" and "Barcode 2".')
parser.add_argument('--out_file', type=str, default='spatial_location_barcodes.csv', help='Outfile to write the results into.')
parser.add_argument('--n_clip_start_1', type=int, default=15, help='Number of characters to clip from the start of each barcode 1 sequence.')
parser.add_argument('--n_clip_end_1', type=int, default=15, help='Number of characters to clip from the end of each barcode 1 sequence.')
parser.add_argument('--n_clip_start_2', type=int, default=24, help='Number of characters to clip from the start of each barcode 2 sequence.')
parser.add_argument('--n_clip_end_2', type=int, default=15, help='Number of characters to clip from the end of each barcode 2 sequence.')
parser.add_argument('--location_id_prefix', type=str, default='A_B', help='Location prefixes separated by an underscore.')

def reverse_complement(sequence):
    # Define the complement mapping
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Generate the complement sequence
    complement_sequence = ''.join(complement[base] for base in sequence)

    # Reverse the complement sequence
    rc_sequence = complement_sequence[::-1]
    
    return rc_sequence

args = parser.parse_args()

# Change directory
os.chdir('/lustre/scratch/kiviaho/dbitseq_atac')

# Assign variables from arguments
xlsx_sheet = args.xlsx_sheet
out_file = args.out_file
n_clip_start_1 = args.n_clip_start_1
n_clip_end_1 = args.n_clip_end_1
n_clip_start_2 = args.n_clip_start_2
n_clip_end_2 = args.n_clip_end_2

location_id_prefix = args.location_id_prefix

id1 = location_id_prefix.split('_')[0]
id2 = location_id_prefix.split('_')[1]

# Read the Excel file
raw_barcodes = pd.read_excel(xlsx_sheet)

# Extract the relevant base sequence and annotate it with a location for barcode 1
barcode1_dict = {}
for i, s in enumerate(raw_barcodes['Barcode 1']):
    barcode1_dict[id1 + str(i + 1)] = s[n_clip_start_1:-n_clip_end_1]

# Extract the relevant base sequence and annotate it with a location for barcode 2
barcode2_dict = {}
for i, s in enumerate(raw_barcodes['Barcode 2']):
    barcode2_dict[id2 + str(i + 1)] = s[n_clip_start_2:-n_clip_end_2]

# Create a list of unique location IDs
unique_pairs = [(a, b) for a in barcode1_dict for b in barcode2_dict]

# Extract the unique sequence combinations by cycling through the oligos
# REMEMBER: It's either RC(barcode 2 + barcode 1) OR RC(barcode 1) + RC(barcode 2). 
# This is due to the ligation direction (backend of seq 2). The first option is implemented here
location_identifier_sequences = [barcode2_dict[tup[1]] + barcode1_dict[tup[0]] for tup in unique_pairs]

# Create a reverse complement of the unique identifiers to match the sequencing libraries
reverse_complement_sequences = [reverse_complement(sequence) for sequence in location_identifier_sequences]

# Create a DataFrame and save it
df = pd.DataFrame({'Sample_ID': ['_'.join(tup) for tup in unique_pairs], 'index': reverse_complement_sequences})
df.to_csv(out_file,index=False)
