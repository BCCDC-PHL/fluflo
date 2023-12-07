#!/usr/bin/env python3 
import os, sys
import re
import argparse
from glob import glob 


def init_parser():
	parser = argparse.ArgumentParser(description="Concatenate sequences from files.")
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument("-s", "--segment", default='NONE', type=str, help="Input folder containing one consensus.fa per isolate, with up to 8 segments per file.")
	group.add_argument("-c", "--concat", action='store_true', help="Input folder containing one consensus.fa per isolate, with up to 8 segments per file.")
	
	parser.add_argument("-i", "--inpath",  help="Input folder containing one consensus.fa per isolate, with up to 8 segments per file.")
	parser.add_argument('-o', "--outpath", type=str, required=True, help="Output file name")
	
	return parser

def parse_fasta(filepath):	
	seqlist = []
	seq = ''
	header = ''
	with open(filepath, 'r') as handle:
		for line in handle.readlines():
			if line[0] == '>':
				if header != '':
					seqlist.append((header, seq))

				header = line.strip().lstrip('>')
				seq = ''
			else:
				seq += line.strip()
		seqlist.append((header, seq))
	
	return seqlist

def write_fasta(seqs, outpath):
	with open(outpath, 'w') as outfile:
		for header, seq in seqs:
			outfile.write(">" + header + '\n')
			outfile.write(seq + '\n')

def concatenate_sequences(filepath):
	sequences = parse_fasta(filepath)

	filepath_fields = re.split("(?:_fluviewer)?[\._]consensus", os.path.basename(filepath))

	if len(filepath_fields) < 2:
		print(f"WARN: Could not extract isolate name from file path {filepath}")
		return None 

	header = filepath_fields[0]

	if len(sequences) != 8:
		print(f"WARN: Skipping {filepath}. Segment sequences are missing")
		return None
	
	concatenated_seq = ''
	for _ , seq in sequences:
		concatenated_seq += seq

	return (header, concatenated_seq)


def extract_sequences(filepath, target_segment):
	sequences = parse_fasta(filepath)

	search_segments = [x for x in sequences if re.search(target_segment, x[0])] 

	if len(search_segments) != 1: 
		print("ERROR: Could not locate a singular sequences that matches the segment of interest. You may need to specify a more advanced regex with --segment")
		print(search_segments)
		return None 

	return search_segments[0]


def main(args):

	if os.path.isfile(args.inpath):
		print("ERROR: Expected directory of consensus.fa files as input.", file=sys.stderr)
		exit(1)
	
	filepaths = glob(os.path.join(args.inpath, '*consensus*'))
	output_sequences = []

	for filepath in filepaths:


		if args.concat:
			seq = concatenate_sequences(filepath)
		else:
			seq = extract_sequences(filepath, args.segment)

		if seq:
			output_sequences.append(seq)

	write_fasta(output_sequences, args.outpath)

if __name__ == "__main__":
	args = init_parser().parse_args()
	main(args)