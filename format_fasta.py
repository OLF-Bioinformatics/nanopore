import sys
import textwrap

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

fasta_dict = dict()

with open(input_fasta, 'r') as in_handle:
    header = ''
    seq = list()
    for line in in_handle:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            if seq:
                fasta_dict[header] = ''.join(seq)  # Store in dictionary
                seq = list()  # empty seq
            header = line
        else:
            seq.append(line)
        # last entry
        fasta_dict[header] = ''.join(seq)

with open(output_fasta, 'w') as out_handle:
    for title, seq in fasta_dict.items():
        out_handle.write('>{}\n{}\n'.format(title, '\n'.join(textwrap.wrap(seq, 80, break_long_words=True))))
