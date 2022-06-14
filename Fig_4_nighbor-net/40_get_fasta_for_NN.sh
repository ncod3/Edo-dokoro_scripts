# Retain biallelic SNPs without missings
#
# STEP 1 : Retain biallelic SNPs                  (bcftools view -m 2 -M 2 -v snps)
# STEP 2 : Remove the positions with missing      (bcftools view -i 'F_MISSING=0')
#
# vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.vcf.gz : original VCF containing all marker information

bcftools view -m 2 \
              -M 2 \
              -v snps \
              vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.vcf.gz | \
bcftools view -i 'F_MISSING=0' \
              -O z \
> vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_missing.m2M2vsnps.vcf.gz


# Convert VCF to genotype matrix
#
# The Rust program 'vcf2gt' is in the directory '/Fig_4_nighbor-net/vcf2gt/target/release'
# Since Rust is compile programming language, you may need to recompile the program in your environment.
# 'vcf2gt' convert '0/0' -> '0', '0/1' -> '1', '1/1' -> 2, and './.' -> '9'.
#
# Each column of vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_missing.m2M2vsnps.vcf.gz.gt indicates: 
# 1.  Chromosome name
# 2.  Position
# 3.  Genotype of dt_kita1  (D. tokoro)
# 4.  Genotype of dt_waka1  (D. tokoro)
# 5.  Genotype of dt_mzawa  (D. tokoro)
# 6.  Genotype of ed_hachi1 (Edo-dokoro)
# 7.  Genotype of tn_yuga   (D. tenuipes)
# 8.  Genotype of tn_utsu   (D. tenuipes)
# 9.  Genotype of tn_mie1   (D. tenuipes)
# 10. Genotype of dq_mie    (D. quinqueloba)

vcf2gt vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_missing.m2M2vsnps.vcf.gz \
     > vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_missing.m2M2vsnps.vcf.gz.gt


# Convert genotype matrix to FASTA
#
# The homemade python script 'vcf2fa' depends on pandas library.
# [usage] vcf2fa.py <vcf.gz> (<No. acceptable missing>)
# In this analysis, we didn't allow any missings.

vcf2fa.py vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_missing.m2M2vsnps.vcf.gz.gt \
        > vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_missing.m2M2vsnps.vcf.gz.fasta

