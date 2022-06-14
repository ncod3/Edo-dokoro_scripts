
# Get genotype table for chloroplast
# Output file 'genotype_table_CP.txt' is in '61_genotype_tables'
# Each column of 'genotype_table_CP.txt' indicates:
# 1.  Contig name (chloroplast)
# 2.  Position
# 3.  REF (nucleotide of reference genome)
# 4.  ALT (nucleotide of mutation)
# 5.  Genotype of dt_kita1   (D. tokoro)
# 6.  Genotype of dt_waka1   (D. tokoro)
# 7.  Genotype of dt_mzawa   (D. tokoro)
# 8.  Genotype of ed_hachi1  (Edo-dokoro)
# 9.  Genotype of tn_yuga    (D. tenuipes)
# 10. Genotype of tn_utsu    (D. tenuipes)
# 11. Genotype of tn_mie1    (D. tenuipes)
# 12. Genotype of dq_mie     (D. quinqueloba)

bcftools view -m 2 \
              -M 2 \
              -i 'F_MISSING=0' \
              vcf_pp_CPMT_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.vcf.gz | \
bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' | \
grep -v '0/1' | \
grep 'CP' > genotype_table_CP.txt



# Get genotype table for mitochondria
# Output file 'genotype_table_MT.txt' is in '61_genotype_tables'
# Each column of 'genotype_table_MT.txt' indicates:
# 1.  Contig name (mitochondria)
# 2.  Position
# 3.  REF (nucleotide of reference genome)
# 4.  ALT (nucleotide of mutation)
# 5.  Genotype of dt_kita1   (D. tokoro)
# 6.  Genotype of dt_waka1   (D. tokoro)
# 7.  Genotype of dt_mzawa   (D. tokoro)
# 8.  Genotype of ed_hachi1  (Edo-dokoro)
# 9.  Genotype of tn_yuga    (D. tenuipes)
# 10. Genotype of tn_utsu    (D. tenuipes)
# 11. Genotype of tn_mie1    (D. tenuipes)
# 12. Genotype of dq_mie     (D. quinqueloba)

bcftools view -m 2 \
              -M 2 \
              -i 'F_MISSING=0' \
              vcf_pp_CPMT_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.vcf.gz | \
bcftools query -f '%CHROM %POS %REF %ALT [%GT ]\n' | \
grep -v '0/1' | \
grep 'MT' > genotype_table_MT.txt



# Convert genotype table to FASTA format
# The python script 'get_CPMT_fasta.py' and output files are in '62_fasta_files'

get_CPMT_fasta.py genotype_table_CP.txt > SNPs_CP.fasta
get_CPMT_fasta.py genotype_table_MT.txt > SNPs_MT.fasta



# Reconstruct organellar trees
# Output files are in '63_trees/CP' and '63_trees/MT'
# Run IQ-TREE for chloroplast markers

iqtree -s SNPs_CP.fasta \
       -m MFP \
       -st DNA \
       -bb 1000

# Run IQ-TREE for mitochondria markers

iqtree -s SNPs_MT.fasta \
       -m MFP \
       -st DNA \
       -bb 1000
