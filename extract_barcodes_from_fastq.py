import argparse
from Bio import SeqIO
import gzip
import pandas as pd
from tqdm import tqdm
import subprocess


index1_start = 15  # 15bp spacer
index1_end = index1_start + 10  # 10bp spatial barcode 1
index2_start = index1_end + 30  # 30bp spacer 
index2_end = index2_start + 10  # 10bp spatial barcode 2

def extract_spatial_barcode(seq_record):
    seq_record1 = seq_record[index1_start:index1_end]
    seq_record2 = seq_record[index2_start:index2_end]
    return seq_record1 + seq_record2

def count_rows_in_fastq_gz(file_path):
    with gzip.open(file_path, 'rt') as file:
        row_count = sum(1 for _ in file)
    return int(row_count / 4)  # 4 rows per sequence in a fastq file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Trim a fastq.gz file containing the spatial barcode readout.")
    parser.add_argument('--data_dir', required=True, help='Directory containing the input file')
    parser.add_argument('--input_file', required=True, help='Name of the input fastq.gz file')
    parser.add_argument('--output_file', required=True, help='Name of the output fastq.gz file')
    parser.add_argument('--valid_barcodes_file', required=True, help='CSV file containing valid barcodes. Barcodes should be in "index" column.')

    args = parser.parse_args()

    data_dir = args.data_dir
    input_file = args.input_file
    output_file = args.output_file
    valid_barcodes_file = args.valid_barcodes_file

    barcode_reference = pd.read_csv(valid_barcodes_file)
    valid_barcodes_set = set(barcode_reference['index'].astype(str))  # Convert to set for faster lookup
    n_total = count_rows_in_fastq_gz(data_dir + input_file)
    n_correct = 0
    progress = 0

    with gzip.open(data_dir + input_file, "rt") as index_file, gzip.open(data_dir + output_file, "wt") as output_handle:
        for record in tqdm(SeqIO.parse(index_file, "fastq"), total=n_total):
            # Trim the reads
            record = extract_spatial_barcode(record)

            if str(record.seq) in valid_barcodes_set:
                n_correct += 1

            progress += 1

            if progress % int(1e6) == 0:
                # Print out the match rate every million reads
                print('1 to 1 barcode match rate: {:.2%} ({:d} million reads processed)'.format(n_correct / progress, progress // int(1e6)))

            # Write the trimmed record to the output file
            SeqIO.write(record, output_handle, "fastq")

    # Remove the original file so that it doesn't clash when processing with cell ranger
    subprocess.run("rm " + data_dir + input_file, shell=True)
