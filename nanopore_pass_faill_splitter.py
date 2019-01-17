#!/usr/local/env python3.6

__version__ = '0.3'
__author__ = ['duceppemo', 'chmaraj']

import gzip
# from math import log
import os
import pathlib
from multiprocessing import Pool
import functions


class Splitter(object):

    def __init__(self, args):
        self.input_fastq = args.input
        self.output_folder = args.output
        self.threshold = args.qscore
        self.threads = args.threads

        # Set the output folder paths for pass and fail reads
        self.pass_dir = '/'.join([self.output_folder, 'pass'])
        self.fail_dir = '/'.join([self.output_folder, 'fail'])

        # Create the folders if they don't exist
        pathlib.Path(self.pass_dir).mkdir(parents=True, exist_ok=True)
        pathlib.Path(self.fail_dir).mkdir(parents=True, exist_ok=True)

        # Run
        self.run()

    def run(self):
        fastq_list = []

        if os.path.isfile(self.input_fastq):
            self.is_fastq(self.input_fastq)  # check if right file extension
            fastq_list.append(self.input_fastq)  # just put that file in the list
        elif os.path.isdir(self.input_fastq):
            for root, directories, filenames in os.walk(self.input_fastq):  # look recursively
                for filename in filenames:
                    absolute_path = os.path.join(root, filename)
                    if os.path.isfile(absolute_path) and filename.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                        fastq_list.append(absolute_path)
        else:
            raise Exception("Something is wrong with your input.")

        # check if input_fastq_list is not empty
        if not fastq_list:
            raise Exception("Input fastq file not found!" if os.path.isfile(self.input_fastq)
                            else "No fastq files found in input folder")

        io_tuple_list = []
        for fastq in fastq_list:
            sample_name = os.path.basename(fastq).split(".")[0]  # everything before the fist dot (".")

            # Create the path for the output file
            output_pass_file = ''.join([self.pass_dir, '/', sample_name, '_pass.fastq.gz'])
            output_fail_file = ''.join([self.fail_dir, '/', sample_name, '_fail.fastq.gz'])

            io_tuple_list.append((fastq, output_pass_file, output_fail_file))

        # Run samples in parallel
        # Initiate the pool
        pool = Pool(self.threads)

        # Parse the fastq files in parallel asynchronously
        pool.map_async(self.parse_fastq_parallel, io_tuple_list)

        # Close the pool
        pool.close()
        pool.join()

    def is_fastq(self, fastq):
        if not 'fastq' or 'fq' in fastq:
            raise Exception('Only files with the following extensions are accepted:\
            ".fastq", ".fq", ".fastq.gz" or ".fq.gz" ')

    def close_handle(self, handle):
        handle.close()

    def parse_fastq_parallel(self, io_tuple):
        # Create output file handles
        in_fastq, pass_fastq, fail_fastq = io_tuple

        pass_file_handle = gzip.open(pass_fastq, 'w')
        fail_file_handle = gzip.open(fail_fastq, 'w')

        handles_tuple = (pass_file_handle, fail_file_handle)
        with gzip.open(in_fastq, 'rb', 1024 * 1024) if \
                in_fastq.endswith('gz') else open(in_fastq, 'rb', 1024 * 1024) as f:
            lines = []
            for line in f:
                if not line:  # end of file?
                    continue
                line = line.rstrip()

                # if lines and b'runid' in line:  # new entry -> specific to nanopore
                if len(lines) == 4:
                    self.process_fastq_entry(lines, handles_tuple)
                    lines = []
                lines.append(line)
            self.process_fastq_entry(lines, handles_tuple)

        self.close_handle(pass_file_handle)
        self.close_handle(fail_file_handle)

    def process_fastq_entry(self, my_list, output_list):
        header, seq, extra, qual = my_list  # get each component of list in a variable

        # Average phred score
        phred_list = [letter - 33 for letter in qual]
        # average_phred = -10 * log(sum([10 ** (q / -10) for q in phred_list]) / len(phred_list), 10)
        average_phred = functions.compute_average_quality(phred_list, len(seq))  # cython

        flag = 'pass'
        if average_phred < self.threshold:
            flag = 'fail'

        self.write_n_compress(flag, header, seq, extra, qual, output_list)

    def write_n_compress(self, flag, header, seq, extra, qual, output_list):
        flag_dict = {'pass': output_list[0],
                     'fail': output_list[1]}

        flag_dict[flag].write(header + b'\n' + seq + b'\n' + extra + b'\n' + qual + b'\n')


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Split a fastq file into pass and fail reads based on average phred score ')
    parser.add_argument('-i', '--input', metavar='sample.fastq.gz|/my/folder',
                        required=True,
                        help='An indidual fastq file or a folder with fastq file(s), gzipped or not')
    parser.add_argument('-o', '--output', metavar='/output/folder/',
                        required=True,
                        help='Output folder'
                             'Output fastq files will be gzipped')
    parser.add_argument('-q', '--qscore', metavar='7', type=float, default=7,
                        required=False,
                        help='Qscore (average phred score) used to split pass reads from fail reads'
                             'Default is 7')
    parser.add_argument('-t', '--threads', metavar='1', type=int, default=1,
                        required=False,
                        help='Threads used in case of multiprocessing files in a folder'
                             'Default is 1')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Splitter(arguments)
