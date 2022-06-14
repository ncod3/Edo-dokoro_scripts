#! /usr/bin/env python3

import sys

# Input genotype table
input_file = sys.argv[1]


with open(input_file) as f:

	seqs = {}
	
	for i in range(0, 8):
		seqs[i] = ''

	for line in f:

		line = line.rstrip('\n')
		cols = line.split(' ')

		ref = cols[2]
		alt = cols[3]

		gts = cols[4:]

		# Extract only homozygous markers from VCF

		for i, gt in enumerate(gts):

			if gt =='0/0':
				
				seqs[i] = seqs[i] + ref

			elif gt == '1/1':

				seqs[i] = seqs[i] + alt

# Output FASTA format
for k, v in seqs.items():

	print('>{}'.format(k))
	print('{}'.format(v))



