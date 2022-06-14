
# STEP 1 : remove the column of dq_mie (D. quinqueloba)  (bcftools view -s ^dq_mie)
# STEP 2 : Retain biallelic SNPs                         (bcftools view -m 2 -M 2 -v snp)
# 
# vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.vcf.gz : original VCF containing all marker information

bcftools view -s ^dq_mie \
              vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.vcf.gz | \
bcftools view -m 2 \
              -M 2 \
              -v snps \
              -O z \
> vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.m2M2vsnps.vcf.gz




# Convert VCF to genotype matrix
# 'vcf2gt' convert '0/0' -> '0', '0/1' -> '1', '1/1' -> 2, and './.' -> '9'
# 'vcf2gt' is in the directory '/Fig_4_nighbor-net/vcf2gt/target/release'
#
# Each column of 'vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.m2M2vsnps.vcf.gz.gt' indicates: 
# 1.  Chromosome name
# 2.  Position
# 3.  Genotype of dt_kita1  (D. tokoro)
# 4.  Genotype of dt_waka1  (D. tokoro)
# 5.  Genotype of dt_mzawa  (D. tokoro)
# 6.  Genotype of ed_hachi1 (Edo-dokoro)
# 7.  Genotype of tn_yuga   (D. tenuipes)
# 8.  Genotype of tn_utsu   (D. tenuipes)
# 9.  Genotype of tn_mie1   (D. tenuipes)

vcf2gt vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.m2M2vsnps.vcf.gz \
     > vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.m2M2vsnps.vcf.gz.gt




# Calculate Fst and pi
# Fst was calculated between D. tokoro and D. tenuipes
# Pi was calculated in each population of D. tokoro and D. tenuipes
#
# Each column of 'calculate_fst_pi_result.txt' indicates: 
# 1 : Chromosome name
# 2 : Fst
# 3 : No. SNPs for Fst
# 4 : Pi in D. tokoro                  (before normalized by mapped sites)
# 5 : Pi in D. tenuipes                (before normalized by mapped sites)
# 6 : pi in D. tokoro and D. tenuipes  (before normalized by mapped sites)
# 7 : No. SNPs for pi in D. tokoro
# 8 : No. SNPs for pi in D. tenuipes
# 9 : No. SNPs for pi in D. tokoro and D. tenuipes

calculate_fst_pi.py vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.m2M2vsnps.vcf.gz.gt \
                  > calculate_fst_pi_result.txt