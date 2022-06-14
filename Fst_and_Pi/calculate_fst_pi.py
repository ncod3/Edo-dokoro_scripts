#! /usr/bin/env python3

import sys


input_file = sys.argv[1]



def count_gts(gts):

    n0 = 2*gts.count(0) + gts.count(1)
    n2 = 2*gts.count(2) + gts.count(1)

    return n0, n2


def calculate_pi(n0, n2):

    n = n0 + n2
    pi = n0*n2/(n*(n - 1)/2)
    pi = pi*n/(n - 1)

    return pi




with open(input_file) as f:

    count_gt = {}
    count_gt['all'] = [0, 0, 0, 0, 0, 0, 0, 0]

    # 0 : Fst
    # 1 : No. SNPs for Fst
    # 2 : Pi in D. tokoro
    # 3 : Pi in D. tenuipes
    # 4 : pi in D. tokoro and D. tenuipes
    # 5 : No. SNPs for Pi in D. tokoro
    # 6 : No. SNPs for Pi in D. tenuipes
    # 7 : No. SNPs for pi in D. tokoro and D. tenuipes

    n_tokoro = 6
    n_hime = 6
    n_all = n_tokoro + n_hime

    for line in f:

        line = line.rstrip('\n')
        cols = line.split(' ')

        chrom_name = cols[0]

        if chrom_name not in count_gt:

            count_gt[chrom_name] = [0, 0, 0, 0, 0, 0, 0, 0]

        gts = cols[2:]

        tokoro_gts = [int(gt) for gt in gts[0:3]]
        edo_gt = [int(gt) for gt in gts[3:4]]
        hime_gts = [int(gt) for gt in gts[4:7]]

        tokoro_miss = tokoro_gts.count(9)
        hime_miss = hime_gts.count(9)

        if tokoro_miss == 0:

            tokoro_0, tokoro_2 = count_gts(tokoro_gts)
            pi = calculate_pi(tokoro_0, tokoro_2)

            if pi != 0:

                count_gt['all'][2] += pi
                count_gt[chrom_name][2] += pi

                count_gt['all'][5] += 1
                count_gt[chrom_name][5] += 1

        
        if hime_miss == 0:

            hime_0, hime_2 = count_gts(hime_gts)
            pi = calculate_pi(hime_0, hime_2)

            if pi != 0:

                count_gt['all'][3] += pi
                count_gt[chrom_name][3] += pi

                count_gt['all'][6] += 1
                count_gt[chrom_name][6] += 1


        if tokoro_miss == 0 and hime_miss == 0:

            tokoro_0, tokoro_2 = count_gts(tokoro_gts)
            hime_0, hime_2 = count_gts(hime_gts)

            all_0 = tokoro_0 + hime_0
            all_2 = tokoro_2 + hime_2

            tokoro_f0 = tokoro_0/(tokoro_0 + tokoro_2)
            hime_f0 = hime_0/(hime_0 + hime_2)

            ave_f0 = (tokoro_f0 + hime_f0)/2

            Ht = 2*ave_f0*(1 - ave_f0)
            Hs = tokoro_f0*(1 - tokoro_f0) + hime_f0*(1 - hime_f0)

            pi = calculate_pi(all_0, all_2)

            if pi != 0:

                count_gt['all'][4] += pi
                count_gt[chrom_name][4] += pi

                count_gt['all'][7] += 1
                count_gt[chrom_name][7] += 1


            if Ht !=0 :

                Fst = (Ht - Hs)/Ht

                count_gt['all'][0] += Fst
                count_gt['all'][1] += 1
                count_gt[chrom_name][0] += Fst
                count_gt[chrom_name][1] += 1


		

for k, v in count_gt.items():

    if v[1] != 0 :

        print(k, v[0]/v[1], v[1], v[2], v[3], v[4], v[5], v[6], v[7], sep='\t')

    else:

        print(k, 0, v[1], v[2], v[3], v[4], v[5], v[6], v[7], sep='\t')
