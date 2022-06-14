#! /usr/bin/env python3

import re
import os
import sys
import gzip
import pandas as pd

print("[usage] vcf2fa.py <vcf.gz> (<No. acceptable missing>)", file=sys.stderr)


if len(sys.argv) < 2:
    print("Please input VCF!", file=sys.stderr)
    sys.exit(1)

vcf_name = sys.argv[1]

if len(sys.argv) < 3:
    N_missing = 0
elif len(sys.argv) == 3:
    N_missing = int(sys.argv[2])
else:
    print("Too many arguments!", file=sys.stderr)
    sys.exit(1)

def iupac(ref, alt, gt):
    if gt == "0/0":
        return ref
    elif gt == "1/1":
        return alt
    elif gt == "0/1":
        if ref == "A" or alt == "A":
            if ref == "T" or alt == "T":
                return "W"
            elif ref == "C" or alt == "C":
                return "M"
            elif ref == "G" or alt == "G":
                return "R"
        elif ref == "C" or alt == "C":
            if ref == "G" or alt == "G":
                return "S"
            elif ref == "T" or alt == "T":
                return "Y"
        else:
            return "K"
    elif gt == "./.":
        return "N"
    else:
        print("{} coudn't be recognized!".format(gt), file=sys.stderr)
        sys.exit(1)

comment = re.compile("[^#]")
with gzip.open(vcf_name, "rt") as vcf:
    with open(vcf_name + ".temp", "w") as temp:
        for line in vcf:
            if comment.match(line):
                line = line.rstrip("\n")
                cols = line.split("\t")
                ref = cols[3]
                alt = cols[4]
                gts = [iupac(ref, alt, field.split(":")[0]) for field in cols[9:]]
                if gts.count("N") <= N_missing:
                    temp.write(" ".join(gts) + "\n")
            else:
                line = line.rstrip("\n")
                cols = line.split("\t")
                sample_names = [os.path.basename(name) for name in cols[9:]]

temp = pd.read_csv(vcf_name + ".temp", sep=" ", header=None)
temp = temp.T

for sample_name, (i, row) in zip(sample_names, temp.iterrows()):
    print(">{}".format(sample_name))
    seq = "".join(list(row))
    print(seq)

os.remove(vcf_name + ".temp")
