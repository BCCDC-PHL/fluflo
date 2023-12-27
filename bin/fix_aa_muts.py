#!/usr/bin/env python3
import os, json
import argparse


def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--aa_json', required=True,  help='aa_muts.json produced by the Augur translate step')	
	parser.add_argument('-n',  '--nt_json', required=True,  help='nt_muts.json produced by the Augur ancestral step. Needed to calculate the length of the full nucleotide sequence.')	
	parser.add_argument('-o', '--outpath', required=True, help='Fixed amino acid mutation node data in JSON format.')	
	
	return parser


def main(args):

	with open(args.aa_json, 'r') as infile:
		aa_data = json.load(infile)

	with open(args.nt_json, 'r') as infile:
		reference_length = len(json.load(infile)['reference']['nuc'])
	
	if 'nuc' not in aa_data['annotations']:
		key = list(aa_data['annotations'].keys())[0]
		name = aa_data['annotations'][key]['seqid']

		aa_data['annotations']['nuc'] = {
		"end": reference_length,
		"seqid": name,
		"start": 1,
		"strand": "+",
		"type": "source"
		}

	with open(args.outpath, 'w') as outfile:
		json.dump(aa_data, outfile, indent=4)

if __name__ == "__main__":
	args = init_parser().parse_args()
	main(args)