
# STEP 1 : Remove the column of dq_mie (D. quinqueloba)  (bcftools view -s ^dq_mie)
# STEP 2 : Retain biallelic SNPs                         (bcftools view -m 2 -M 2 -v snp)
# STEP 3 : Remove the positions with missing             (bcftools view -i 'F_MISSING=0')
# 
# vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.vcf.gz : original VCF containing all marker information

bcftools view -s ^dq_mie \
              vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.vcf.gz | \
bcftools view -m 2 \
              -M 2 \
              -v snps | \
bcftools view -i 'F_MISSING=0' \
              -O z \
> vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.no_missing.m2M2vsnps.vcf.gz



# Convert VCF to genotype matrix
# 'vcf2gt' convert '0/0' -> '0', '0/1' -> '1', '1/1' -> 2, and './.' -> '9'
# 'vcf2gt' is in the directory '/Fig_4_nighbor-net/vcf2gt/target/release'
#
# Each column of 'vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.no_missing.m2M2vsnps.vcf.gz.gt' indicates: 
# 1.  Chromosome name
# 2.  Position
# 3.  Genotype of dt_kita1  (D. tokoro)
# 4.  Genotype of dt_waka1  (D. tokoro)
# 5.  Genotype of dt_mzawa  (D. tokoro)
# 6.  Genotype of ed_hachi1 (Edo-dokoro)
# 7.  Genotype of tn_yuga   (D. tenuipes)
# 8.  Genotype of tn_utsu   (D. tenuipes)
# 9.  Genotype of tn_mie1   (D. tenuipes)

vcf2gt vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.no_missing.m2M2vsnps.vcf.gz \
     > vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.no_missing.m2M2vsnps.vcf.gz.gt



# Fig. 5A
# Count Edo-dokoro genotypes on the positions that D. tokoro and D. tenuipes genotypes are oppositely fixed.
#
# Each column of 'count_edo_gt_type_result' indicates:
# 1. Chromosome name
# 2. Number of homozygous genotype of D. tokoro
# 3. Number of hetero between D. tokoro nad D. tenuipes
# 4. Number of homozygous genotype of D. tenuipes

./count_edo_gt_type.py vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.no_missing.m2M2vsnps.vcf.gz.gt \
                > count_edo_gt_type_result.txt


# Fig. 5B
# Check Edo-dokoro genotypes on the positions that D. tokoro and D. tenuipes genotypes are oppositely fixed.
# Each column of 'check_edo_gt_type_result.txt' indicates:
# 1. Chromosome name
# 2. Position
# 3. Genotype of Edo-dokoro (0: homozygous genotype of D. tokoro, \
#                            1: hetero between D. tokoro nad D. tenuipes, \
#							 2: homozygous genotype of D. tenuipes)

check_edo_gt_type.py vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_dq_mie.no_missing.m2M2vsnps.vcf.gz.gt \
              > check_edo_gt_type_result.txt
