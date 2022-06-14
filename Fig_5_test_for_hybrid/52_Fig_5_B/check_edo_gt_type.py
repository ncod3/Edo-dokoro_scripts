#! /usr/bin/env python3

import sys


input_file = sys.argv[1]

with open(input_file) as f:

    for line in f:

        line = line.rstrip('\n')
        cols = line.split(' ')

        chrom_name = cols[0]
        posi = cols[1]

        gts = cols[2:]

        tokoro_gts = [int(gt) for gt in gts[0:3]]
        edo_gt = [int(gt) for gt in gts[3:4]]
        hime_gts = [int(gt) for gt in gts[4:7]]

        if sum(tokoro_gts) == len(tokoro_gts)*0 and sum(hime_gts) == len(hime_gts)*2:

            print(chrom_name, posi, edo_gt[0], sep='\t')

        elif sum(tokoro_gts) == len(tokoro_gts)*2 and sum(hime_gts) == len(hime_gts)*0:

            print(chrom_name, posi, 2 - edo_gt[0], sep='\t')
