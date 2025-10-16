#!/usr/bin/env python
from Bio import SeqIO
import pandas as pd
import argparse

# Define Function to read diamond oufmt 6 output
def read_diamond(file_path):
    try:
        df = pd.read_csv(file_path, sep="\t", header=None,
                             names=["qseqid", "sseqid", "pident", "length", "mismatch", 
                                    "gapopen", "qstart", "qend", "sstart", "send", 
                                    "evalue", "bitscore"])
        return(df)
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return pd.DataFrame()
    except pd.errors.EmptyDataError:
        print(f"Error: File '{file_path}' is empty.")
        return pd.DataFrame()
    except Exception as e:
        print(f"Error: Failed to parse file '{file_path}'. {str(e)}")
        return pd.DataFrame()


# Define Function to read nucleotide sequence file
def read_fasta(file_path):
    try:
        return SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return []
    except Exception as e:
        print(f"Error: Failed to parse file '{file_path}'. {str(e)}")
        return []
    
def main():
    parser = argparse.ArgumentParser(
                    prog='Extract_Fastas_Diamond.py',
                    description='Python script to extract aligned sequences from the query file of a Diamond outfmt 6 output',
                    epilog='Developed by MsC. Matheus Cosentino, 29.05.2025')

    parser.add_argument('--diamond', help='Diamond output file to map our fasta extraction')
    parser.add_argument('--fasta', help='Query file used in the previously blast file')
    parser.add_argument('--output', help='File with extracted region from alignment')
    args = parser.parse_args() 

    #Read Diamond file of interest
    df = read_diamond(args.diamond)
    if df.empty:
        print("No data loaded from DIAMOND output. Exiting.")
    else:
        print("Diamond File Readed")  # Print only first few rows
        
    #Read the Fasta File of Interest    
    fasta = read_fasta(args.fasta)

    if not fasta:  # Empty list means no sequences
        print("No fasta data loaded. Exiting")
    else:
        print("Fasta File Readed") 

    with open(args.output, "w") as out_fasta:
        for idx, row in df.iterrows():
            qid = row['qseqid']
            start = int(row['qstart'])
            end = int(row['qend'])
            seq_record = fasta[qid]
            if start > end:
                start, end = end, start
                subseq = seq_record.seq[start-1:end].reverse_complement()
                strand = "-"
            else:
                subseq = seq_record.seq[start-1:end]
                strand = "+"
            out_fasta.write(f">{qid}_region_{start}_{end}_{strand}\n{subseq}\n")
            print(f"Extraction complete. Output written to {args.output}")


if __name__ == "__main__":
    main()
