#!/usr/bin/env python3 
import os, sys
from Bio import SeqIO
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

def concatenate_sequences(filepath):
	sequences = list(SeqIO.parse(filepath, "fasta"))
	
	if len(sequences) == 8:
		concatenated_seq = sequences[0]
		for seq in sequences[1:]:
			concatenated_seq += seq

		# rename sequences
		concatenated_seq.id = re.split("_fluviewer|_segment", sequences[0].id)[0]
		concatenated_seq.name = concatenated_seq.id
		concatenated_seq.description = concatenated_seq.id
		return concatenated_seq
	
	print(f"WARN: Skipping {filepath}. Segment sequences are missing")

	return None

def extract_sequences(filepath, target_segment):
	sequences = list(SeqIO.parse(filepath, "fasta"))

	for seq in sequences:
		segment_search = re.search("(?:_fluviewer\||_segment\d_)([A-Z0-9]+)", seq.id)

		if segment_search and segment_search.group(1) == target_segment:
			return seq

	print(f"WARN: Skipping {filepath}. Could not find segment of interest.")
	return None

def main(args):

	if os.path.isfile(args.inpath):
		print("ERROR: Expected directory of consensus.fa files as input.", file=sys.stderr)
		exit(1)
	
	filepaths = glob(os.path.join(args.inpath, '*consensus*'))

	if args.concat:
		output_sequences = []
		for filepath in filepaths:
			concatenated_seq = concatenate_sequences(filepath)
			if concatenated_seq:
				output_sequences.append(concatenated_seq)


	else:
		output_sequences = []
		for filepath in filepaths:
			extracted_seq = extract_sequences(filepath, args.segment)
			if extracted_seq:
				output_sequences.append(extracted_seq)

	SeqIO.write(output_sequences, args.outpath, "fasta")

if __name__ == "__main__":
	args = init_parser().parse_args()
	main(args)