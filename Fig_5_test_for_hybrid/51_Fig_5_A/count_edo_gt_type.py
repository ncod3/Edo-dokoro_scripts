#! /usr/bin/env python3

import sys


input_file = sys.argv[1]

with open(input_file) as f:

    # count_gt['chromosome'] = [N_tokoro_homo_type_SNP, \
    #                           N_hetero_SNP, \
    #                           N_tenuipes_homo_type_SNP]

    count_gt = {}
    count_gt['all'] = [0, 0, 0]


    for line in f:

        line = line.rstrip('\n')
        cols = line.split(' ')

        chrom_name = cols[0]

        if chrom_name not in count_gt:

            count_gt[chrom_name] = [0, 0, 0]

        gts = cols[2:]

        tokoro_gts = [int(gt) for gt in gts[0:3]]
        edo_gt = [int(gt) for gt in gts[3:4]]
        hime_gts = [int(gt) for gt in gts[4:7]]

        if sum(tokoro_gts) == len(tokoro_gts)*0 and sum(hime_gts) == len(hime_gts)*2:

            count_gt['all'][0] += edo_gt.count(0)
            count_gt['all'][1] += edo_gt.count(1)
            count_gt['all'][2] += edo_gt.count(2)

            count_gt[chrom_name][0] += edo_gt.count(0)
            count_gt[chrom_name][1] += edo_gt.count(1)
            count_gt[chrom_name][2] += edo_gt.count(2)


        elif sum(tokoro_gts) == len(tokoro_gts)*2 and sum(hime_gts) == len(hime_gts)*0:

            count_gt['all'][0] += edo_gt.count(2)
            count_gt['all'][1] += edo_gt.count(1)
            count_gt['all'][2] += edo_gt.count(0)

            count_gt[chrom_name][0] += edo_gt.count(2)
            count_gt[chrom_name][1] += edo_gt.count(1)
            count_gt[chrom_name][2] += edo_gt.count(0)
		

for k, v in count_gt.items():

    print(k, v[0], v[1], v[2], sep='\t')
