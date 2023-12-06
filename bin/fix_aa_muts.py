#!/usr/bin/env python3
import os, json
import argparse


def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument( 'json',  help='aa_muts.json produced by the Augur translate step')	
	parser.add_argument( 'outpath',  help='Fixed amino acid mutation node data in JSON format.')	
	
	return parser


def main(args):

	with open(args.json, 'r') as infile:
		data = json.load(infile)


	name = data['annotations']['HA']['seqid']

	data['annotations']['nuc'] = {
      "end": 13590,
      "seqid": name,
      "start": 1,
      "strand": "+",
      "type": "source"
    }

	with open(args.outpath, 'w') as outfile:
		json.dump(data, outfile, indent=4)

if __name__ == "__main__":
	args = init_parser().parse_args()
	main(args)