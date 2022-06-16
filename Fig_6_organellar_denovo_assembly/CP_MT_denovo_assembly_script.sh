#! /bin/sh

threads=36
genome_size=400m

#0.short_fastq
#|-- qc_1.fastq.gz
#`qc_2.fastq.gz

#1.ONT_fastq
#`-- ONT.fastq.gz

#necat: de novo assembly
mkdir 2.assemble
cd 2.assemble

necat bridge config.txt

#racon: polish with nanopore long read
mkdir ../3.racon
cd ../3.racon

minimap2 -t ${threads} -x map-ont --secondary=no \
    ../2.assemble/sample/6-bridge_contigs/polished_contigs.fasta \
    ../1.ONT_fastq/ONT.fastq.gz \
    > long_reads_mapped.paf
racon -t ${threads} \
    ../1.ONT_fastq/ONT.fastq.gz \
    long_reads_mapped.paf \
    ../2.assemble/sample/6-bridge_contigs/polished_contigs.fasta \
    > necat_racon.fasta

#medaka: correct misassembly
mkdir ../4.medaka
cd ../4.medaka

medaka_consensus \
    -i ../1.ONT_fastq/ONT.fastq.gz \
    -d ../3.racon/necat_racon.fasta \
    -o racon_medaka \
    -t ${threads} \
    -m r941_prom_hac_g507

#hypo: polish with Illumina short reads
mmkdir ../5.hypo
cd ../5.hypo

bwa-mem2 index ../4.medaka/racon_medaka/consensus.fasta

bwa-mem2 mem -M -a -t ${threads} \
    ../4.medaka/racon_medaka/consensus.fasta \
    ../0.short_fastq/qc_1.fastq.gz ../0.short_fastq/qc_2.fastq.gz \
    | samtools sort -@ ${threads} -m 20G -o short_1.bam && samtools index -@ ${threads} short_1.bam

coverm contig -t ${threads} -b short_1.bam > coverm_1_contig.txt

seqkit fx2tab -ln ../4.medaka/racon_medaka/consensus.fasta | sort -k 2 -nr > sort.medaka.fasta.txt
seqkit fx2tab -ln ../4.medaka/racon_medaka/consensus.fasta \
    | sort -k 2 -nr | head -n 10 | cut -f 1 \
    > sort.medaka.fasta.top10contig.txt

top10contig=$(cat sort.medaka.fasta.top10contig.txt | tr '\n' ' ')
for i in ${top10contig}
do
  cat coverm_1_contig.txt | grep "$i" >> top10contig.cov1.txt
done
COV1=$(cat top10contig.cov1.txt | awk '{sum+=$2} END {print sum/NR}')
rm -rf sort.medaka.fasta.top10contig.txt

hypo -d ../4.medaka/racon_medaka/consensus.fasta -r @names.txt -s ${genome_size} -c ${COV1} -b short_1.bam -p 12 -t ${threads} -o hypo_1.fasta

# 2nd hypo

bwa-mem2 index hypo_1.fasta

bwa-mem2 mem -M -a -t ${threads} \
    hypo_1.fasta \
    ../0.short_fastq/qc_1.fastq.gz ../0.short_fastq/qc_2.fastq.gz \
    | samtools sort -@ ${threads} -m 20G -o hypo1.bam && samtools index -@ ${threads} hypo1.bam

coverm contig -t ${threads} -b hypo1.bam > coverm_2_contig.txt

seqkit fx2tab -ln hypo_1.fasta | sort -k 2 -nr > sort.hypo_1.fasta.txt
seqkit fx2tab -ln hypo_1.fasta \
    | sort -k 2 -nr | head -n 10 | cut -f 1 \
    > sort.medaka.fasta.top10contig.txt

top10contig=$(cat sort.medaka.fasta.top10contig.txt | tr '\n' ' ')
for i in ${top10contig}
do
  cat coverm_1_contig.txt | grep "$i" >> top10contig.cov2.txt
done
COV2=$(cat top10contig.cov2.txt | awk '{sum+=$2} END {print sum/NR}')
rm -rf sort.medaka.fasta.top10contig.txt

hypo -d hypo_1.fasta -r @names.txt -s ${genome_size} -c ${COV2} -b hypo1.bam -p 12 -t ${threads} -o hypo_2.fasta

# purge_haplotigs
mkdir ../6.purge_haplotigs
cd ../6.purge_haplotigs

minimap2 -t 24 -ax map-ont ../5.hypo/hypo_2.fasta \
  ../1.ONT_fastq/*.fastq.gz \
  --secondary=no \
  | samtools sort -O BAM -o hypo_2_aligned.bam -T tmp.ali

seqkit fx2tab -lni hypo_2.fasta | sort > length.tsv
coverm contig -t 8 -b hypo_2_aligned.bam | sort > depth.tsv
join length.tsv depth.tsv | tr ' ' '\t' > length_depth.tsv

purge_haplotigs hist -b hypo_2_aligned.bam -g ../5.hypo/hypo_2.fasta -t 24
# Check hypo_2_aligned.bam.histogram.png and set parameter
purge_haplotigs cov -i hypo_2_aligned.bam.gencov -l 5 -m 25 -h 10000
purge_haplotigs purge -g ../5.hypo/hypo_2.fasta -c coverage_stats.csv -d -b hypo_2_aligned.bam -t 12
coverm contig -b hypo_2_aligned.bam > coverm.aligned.bam.txt

cat curated.reassignment.tsv| sort > sorted_reassignment.tsv
join length_depth.tsv sorted_reassignment.tsv | column -t > length_depth_reassignment.tsv

cat length_depth_reassignments.tsv | awk '$3>430{print $1}' > deep_contigs.list
seqkit grep -f deep_contigs.list ../5.hypo/hypo_2.fasta > deep_contigs.fasta

# organelle genome assembly
# cp: chloroplast
# mt: mitochondria
mkdir ../7.organelle
cd ../7.organelle

mkdir blast
cd blast

cp ../../6.purge_haplotigs/deep_contigs.fasta .

makeblastdb -in deep_contigs.fasta -dbtype nucl

blastn -query Dioscorea_rotundata_plastid_complete_genome.NC_024170.1_155406.fasta \
-db deep_contigs.fasta \
-outfmt 6 \
-max_target_seqs 10 -max_hsps 1 -perc_identity 98 \
> cp_blast.out

cat cp_blast.out | awk '{print $2}' | sort | uniq > cp_contigs.list

blastn -query Dioscorea_rotundata_mitochondrial_DNA_contig.fasta \
-db deep_contigs.fasta \
-outfmt 6 \
-max_target_seqs 10 -max_hsps 1 -perc_identity 98 \
> mt_blast.out

cat mt_blast.out | awk '{print $2}' | sort | uniq > mt_contigs.list

cd ../

cat ../cp_contigs.list | while read line
do samtools view ../6.purge_haplotigs/hypo_2_aligned.bam $line -bh | samtools bam2fq
done | seqkit seq -m 10000 -M 100000 > cp.fastq

seqkit stats -a cp.fastq
#downsize cp.fastq
seqkit sample -p 0.008 cp.fastq > cp4flye_p0008.fastq
seqkit stats -a cp4flye_p0008.fastq

flye --nano-raw cp4flye.fastq -g 150k -t 12 -o cp_flye

cat ../mt_contigs.list | while read line
do samtools view ../6.purge_haplotigs/hypo_2_aligned.bam $line -bh | samtools bam2fq
done | seqkit seq -m 10000 -M 100000 > mt.fastq

seqkit stats -a mt.fastq
#downsize mt.fastq
seqkit sample -p 0.12 mt.fastq > mt4flye_p012.fastq
seqkit stats -a mt4flye_p012.fastq

flye --nano-raw mt4flye.fastq -g 500k -t 12 -o mt_flye

#re-run the bridge stage of NECAT
cd ../2.assemble/sample
cp -r 4-fsa original_4-fsa
cd original_4-fsa

mkdir 3.racon
cd 3.racon

minimap2 -t ${threads} -x map-ont --secondary=no \
    ../contigs.fasta \
    ../../../../1.ONT_fastq/ONT.fastq.gz \
    > contigs_fasta.paf
racon -t ${threads} \
    ../../../../1.ONT_fastq/ONT.fastq.gz \
    contigs_fasta.paf \
    ../contigs.fasta \
    > contigs_racon.fasta

mkdir ../4.medaka
cd ../4.medaka

medaka_consensus \
    -i ../../../../1.ONT_fastq/ONT.fastq.gz \
    -d ../3.racon/contigs_racon.fasta \
    -o racon_medaka \
    -t 72 \
    -m r941_prom_hac_g507

mkdir ../6.purge_haplotig
cd ../6.purge_haplotig

minimap2 -t ${threads} -ax map-ont --secondary=no \
  ../4.medaka/racon_medaka/consensus.fasta \
  ../../../../1.ONT_fastq/ONT.fastq.gz \
  | samtools sort -@ ${threads} -m 20G -o consensus.bam && samtools index -@ ${threads} consensus.bam

purge_haplotigs hist -b consensus.bam -g ../4.medaka/racon_medaka/consensus.fasta -t ${threads}

seqkit fx2tab -lni ../4.medaka/racon_medaka/consensus.fasta | sort > length.tsv
coverm contig -b consensus.bam -t 24 | sort > depth.tsv
join length.tsv depth.tsv > length_depth.tsv

purge_haplotigs cov -i consensus.bam.gencov -l 10 -m 25 -h 200
purge_haplotigs purge -g ../4.medaka/racon_medaka/consensus.fasta -c coverage_stats.csv -d -b consensus.bam -t 24

cat curated.reassignments.tsv | sort > sorted_reassignments.tsv
join length_depth.tsv sorted_reassignments.tsv | column -t > length_depth_reassignments.tsv
cat length_depth_reassignments.tsv | grep "REPEAT" | awk '$2>20000 && $3>28 && $3<59 {print $1}' > rescue.list
seqkit grep -f rescue.list ../4.medaka/racon_medaka/consensus.fasta > rescue.fasta

cat curated.fasta rescue.fasta \
../../../../7.organelle/cp_flye/assembly.fasta \
../../../../7.organelle/mt_flye/assembly.fasta > for_bridge.fasta

cp for_bridge.fasta ../../4-fsa/contigs.fasta

cd ../../../

necat bridge config.txt

#hypo: polish with Illumina short reads
mkdir ../5.hypo_after_purge_haplotigs_necat_bridge
cd ../5.hypo_after_purge_haplotigs_necat_bridge

cp ../2.assemble/Dt_q10_l1k/6-bridge_contigs/polished_contigs.fasta polished_contigs_CPMT.fasta
#Manual curation were made, such as changing contig names of CP and MT.

bwa-mem2 index polished_contigs_CPMT.fasta

bwa-mem2 mem -M -a -t ${threads} \
    polished_contigs_CPMT.fasta \
    ../0.short_fastq/qc_1.fastq.gz ../0.short_fastq/qc_2.fastq.gz \
    | samtools sort -@ ${threads} -m 20G -o short_1.bam && samtools index -@ ${threads} short_1.bam

coverm contig -t ${threads} -b short_1.bam > coverm_1_contig.txt

seqkit fx2tab -ln polished_contigs_CPMT.fasta | sort -k 2 -nr > sort.medaka.fasta.txt
seqkit fx2tab -ln polished_contigs_CPMT.fasta \
    | sort -k 2 -nr | head -n 10 | cut -f 1 \
    > sort.medaka.fasta.top10contig.txt

top10contig=$(cat sort.medaka.fasta.top10contig.txt | tr '\n' ' ')
for i in ${top10contig}
do
  cat coverm_1_contig.txt | grep "$i" >> top10contig.cov1.txt
done
COV1=$(cat top10contig.cov1.txt | awk '{sum+=$2} END {print sum/NR}')
rm -rf sort.medaka.fasta.top10contig.txt

hypo -d polished_contigs_CPMT.fasta -r @names.txt -s ${genome_size} -c ${COV1} -b short_1.bam -p 12 -t ${threads} -o hypo_1.fasta

# 2nd hypo

bwa-mem2 index hypo_1.fasta

bwa-mem2 mem -M -a -t ${threads} \
    hypo_1.fasta \
    ../0.short_fastq/qc_1.fastq.gz ../0.short_fastq/qc_2.fastq.gz \
    | samtools sort -@ ${threads} -m 20G -o hypo1.bam && samtools index -@ ${threads} hypo1.bam

coverm contig -t ${threads} -b hypo1.bam > coverm_2_contig.txt

seqkit fx2tab -ln hypo_1.fasta | sort -k 2 -nr > sort.hypo_1.fasta.txt
seqkit fx2tab -ln hypo_1.fasta \
    | sort -k 2 -nr | head -n 10 | cut -f 1 \
    > sort.medaka.fasta.top10contig.txt

top10contig=$(cat sort.medaka.fasta.top10contig.txt | tr '\n' ' ')
for i in ${top10contig}
do
  cat coverm_1_contig.txt | grep "$i" >> top10contig.cov2.txt
done
COV2=$(cat top10contig.cov2.txt | awk '{sum+=$2} END {print sum/NR}')
rm -rf sort.medaka.fasta.top10contig.txt

hypo -d hypo_1.fasta -r @names.txt -s ${genome_size} -c ${COV2} -b hypo1.bam -p 12 -t ${threads} -o hypo_2.fasta

#hypo2.fasta is a final product of de novo assembly, including CP and MT sequences. 


#### Reference ####
#### Niato K. Nanopore sequencing for plant genome (in Japanese).
#### in Arakawa K. and Miyamoto M. (eds.),
#### Guide for the long-read WET & DRY analysis (published in Japanese).
#### YODOSYA Co., LTD. Tokyo, Japan. pp.177-199. (2021) ISBN978-4-7581-2253-5

