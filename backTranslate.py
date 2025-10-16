#!/usr/bin/python
#Sript to back translate aligned aminoacid sequences to nucleotide sequences

# Import packages to be used in the script
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Define Function to read aminoacid alignment file
def parse_amino_acid_alignment(file_path):
    try:
        return list(SeqIO.parse(file_path, "fasta"))
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return []
    except Exception as e:
        print(f"Error: Failed to parse file '{file_path}'. {str(e)}")
        return []

# Define Function to read nucleotide sequence file
def parse_nucleotide_sequences(file_path):
    try:
        return list(SeqIO.parse(file_path, "fasta"))
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return []
    except Exception as e:
        print(f"Error: Failed to parse file '{file_path}'. {str(e)}")
        return []

# Define Function of Back Translation
def back_translate(amino_acid_sequences, nucleotide_sequences, output_file):
    with open(output_file, "w") as out_file:
        for amino_acid_sequence in amino_acid_sequences:
            amino_acid_seq = str(amino_acid_sequence.seq)
            for nucleotide_sequence in nucleotide_sequences:
                if nucleotide_sequence.id == amino_acid_sequence.id:
                    nucleotide_seq = str(nucleotide_sequence.seq)
                    chain = ""
                    for aa in amino_acid_seq:
                        if aa != "-":
                            codon = nucleotide_seq[:3].ljust(3, "-")
                            nucleotide_seq = nucleotide_seq[3:]
                            chain += codon
                        else:
                            chain += "---"
                    sequence_object = Seq(chain)
                    record = SeqRecord(sequence_object, id=nucleotide_sequence.id, description="")
                    SeqIO.write(record, out_file, "fasta")

def main():
    parser = argparse.ArgumentParser(
      prog='backTranslate.py',
      description='This script automates the process of back-translating a protein alignment into a corresponding nucleotide alignment',
      epilog='Developed by MsC. Matheus Cosentino, 30.06.2025')
    parser.add_argument("-nt", dest="nucleotide_file", help="nucleotide sequence file", type=str)
    parser.add_argument("-ali", dest="amino_acid_file", help="amino acid multiple alignment file", type=str)
    parser.add_argument("-out", dest="output_file", help="output file for nucleotide multiple alignment", type=str)
    args = parser.parse_args()

    amino_acid_sequences = parse_amino_acid_alignment(args.amino_acid_file)
    nucleotide_sequences = parse_nucleotide_sequences(args.nucleotide_file)

    if amino_acid_sequences and nucleotide_sequences:
        back_translate(amino_acid_sequences, nucleotide_sequences, args.output_file)
    else:
        print("Error: Failed to parse input files.")

if __name__ == "__main__":
    main()
    main()
