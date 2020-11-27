import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
import textwrap

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

fasta_dict = dict()

with open(input_fasta, 'r') as in_handle:
    for title, seq in SimpleFastaParser(in_handle):
        fasta_dict[title] = seq

with open(output_fasta, 'w') as out_handle:
    for title, seq in fasta_dict.items():
        out_handle.write('>{}\n{}\n'.format(title, '\n'.join(textwrap.wrap(seq, 80, break_long_words=True))))
        
