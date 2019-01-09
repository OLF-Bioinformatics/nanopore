#!/usr/local/env python3.6

__version__ = '0.1'
__author__ = ['duceppemo', 'chmaraj']


import gzip
from math import log
import os


class Splitter(object):

    def __init__(self, args):
        self.input_fastq = args.input
        self.output_folder = args.output
        self.threshold = args.qscore

        sample_name = os.path.basename(self.input_fastq).split('.')[0]
        pass_file = ''.join([self.output_folder, '/', sample_name, '_pass.fastq.gz'])
        fail_file = ''.join([self.output_folder, '/', sample_name, '_fail.fastq.gz'])
        self.pass_file_handle = gzip.open(pass_file, 'w')
        self.fail_file_handle = gzip.open(fail_file, 'w')

        self.run()

    def run(self):
        self.parse_fastq(self.input_fastq)
        self.close_handle(self.pass_file_handle)
        self.close_handle(self.fail_file_handle)

    def parse_fastq(self, fastq):

        with gzip.open(fastq, 'rb', 1024 * 1024) if fastq.endswith('gz') else open(fastq, 'rb', 1024 * 1024) as f:
            lines = []
            for line in f:
                if not line:  # end of file?
                    continue
                line = line.rstrip()

                # if lines and b'runid' in line:  # new entry -> specific to nanopore
                if len(lines) == 4:
                    self.process_fastq_entry(lines)
                    lines = []
                lines.append(line)
            self.process_fastq_entry(lines)  # For the last entry

    def process_fastq_entry(self, my_list):
        header, seq, extra, qual = my_list  # get each component of list in a variable

        # Average phred score
        phred_list = [letter - 33 for letter in qual]
        average_phred = -10 * log(sum([10 ** (q / -10) for q in phred_list]) / len(phred_list), 10)

        flag = 'pass'
        if average_phred < 7:
            flag = 'fail'

        self.write_n_compress(flag, header, seq, extra, qual)

    def write_n_compress(self, flag, header, seq, extra, qual):
        flag_dict = {'pass': self.pass_file_handle,
                     'fail': self.fail_file_handle
                     }

        flag_dict[flag].write(header + b'\n' + seq + b'\n' + extra + b'\n' + qual + b'\n')

    def close_handle(self, handle):
        handle.close()


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Split a fastq file into pass and fail reads based on average phred score ')
    parser.add_argument('-i', '--input', metavar='sample.fastq.gz',
                        required=False,
                        help='Input folder with fastq file(s),gzipped or not')
    parser.add_argument('-o', '--output', metavar='/output/folder/',
                        required=True,
                        help='Output folder'
                             'Output fastq files will be gzipped')
    parser.add_argument('-q', '--qscore', metavar='7', type=float, default=7,
                        required=False,
                        help='Qscore (average phred score) used to split pass reads from fail reads'
                             'Default is 7')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Splitter(arguments)
