#! /bin/bash

bams=$@

cpu=10
fasta="DTokoro_PseudoChromosome.vP2.rev01_w1808.fasta"
region="all"
outvcf="vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.vcf.gz"

if [ "x$region" = "xall" ]; then
    region=""
else
    region="-f $region"
fi

bcftools \
    mpileup \
    -a DP,AD,ADF,ADR,SP \
    -B \
    -q 10 \
    -Q 13 \
    -C 50 \
    -I \
    --ignore-RG \
    -O u \
    $region \
    -f $fasta \
		$bams \
\
| bcftools \
	call \
	-vm \
	-f GQ,GP \
	-O u \
\
| bcftools \
	filter \
	-i 'INFO/MQ>=10' \
	-O z \
	--threads $cpu \
	-o $outvcf


