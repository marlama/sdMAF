# pipeline used for sdMAF application in Austism WGS data:

### After applying the Harmonization pipeline from (Leal., et al 2022) we needed to change the CHR code before use the sdMAF algorithm:

### Filter Cases 
plink2 --vcf ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR.vcf.gz --make-bed --out ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR

awk  '{if ($2 == 1) {print "0""\t"$1}}' XWAS_Covariants.txt > XWAS_CasesList.txt

plink --bfile ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR --keep XWAS_CasesList.txt --make-bed --out ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_OnlyCases

418652 variants and 6873 people


### Add Info
awk  '{if ($2 == 1) {print "0""\t"$1"\t"$4}}' XWAS_Covariants.txt > XWAS_CasesSexInfo.txt

awk  '{if ($2 ==1) {print "0""\t"$1"\t"$3"\t"$1}}' XWAS_Covariants.txt > XWAS_CasesPopInfo.txt

awk  '{if ($2 == 0) {print "0""\t"$1"\t"$4}}' XWAS_Covariants.txt > XWAS_ControlsSexInfo.txt

awk  '{if ($2 ==0) {print "0""\t"$1"\t"$3"\t"$1}}' XWAS_Covariants.txt > XWAS_ControlsPopInfo.txt

for pop in Cases Controls ; do plink --bfile ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop} --update-sex XWAS_${pop}SexInfo.txt --make-bed --out ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_sexInfo ; done
Cases 418652 variants and 6873 people
Controls 418652 variants and 8981 people

for pop in Cases Controls ; do plink --bfile ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_sexInfo --update-ids XWAS_${pop}PopInfo.txt --make-bed --out ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_Info ; done
Cases 418652 variants and 6873 people
Controls 418652 variants and 8981 people

### Change chr connotation
for pop in Cases Controls ; do awk '{print$2"\t""22"}' ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_Info_PAR.bim > ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_Info_PAR_UpdateID.txt ; done

for pop in Cases Controls ; do plink --bfile ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_Info --update-chr ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_Info_PAR_UpdateID.txt --make-bed --out ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_Info_UpdatedCHR ; done 

### Generate geno counts
for pop in Cases Controls ; do plink2 --bfile ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_Info_UpdatedCHR  --geno-counts --keep-males --out chrX_MalesGCount_ASD${pop}_XWAS1_CHR22 ; done
Cases
1234 samples removed due to sex filter(s).
5639 samples (0 females, 5639 males; 5639 founders) remaining after main

Controls
5070 samples removed due to sex filter(s).
3911 samples (0 females, 3911 males; 3911 founders) remaining after main

for pop in Cases Controls ; do plink2 --bfile ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_Info_UpdatedCHR --geno-counts --keep-females --out chrX_FemalesGCount_ASD${pop}_XWAS1_CHR22 ; done
Cases
5639 samples removed due to sex filter(s).
1234 samples (1234 females, 0 males; 1234 founders) remaining after main

Controls
3911 samples removed due to sex filter(s).
5070 samples (5070 females, 0 males; 5070 founders) remaining after main

### Formated geno counts
for pop in Cases Controls; do awk '{if (NR==1) {print} if (NR>1) {print "23""\t"$2"\t"$3"\t"$4"\t"$8"\t"$9"\t"$7"\t"$5"\t"$6"\t"$10}}' chrX_MalesGCount_ASD${pop}_XWAS1_CHR22.gcount > chrX_MalesGCount_ASD${pop}_XWAS1_CHR22_Formated1.gcount ; done

for pop in Cases Controls; do awk '{if (NR==1) {print} if (NR>1) {print "23""\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}' chrX_FemalesGCount_ASD${pop}_XWAS1_CHR22.gcount > chrX_FemalesGCount_ASD${pop}_XWAS1_CHR22_Formated1.gcount ; done


### Run sdMAF
module load R/4.2.2_tcltk
for pop in Cases Controls ; do Rscript sdMAF.R --female chrX_FemalesGCount_ASD${pop}_XWAS1_CHR22_Formated1.gcount --male chrX_MalesGCount_ASD${pop}_XWAS1_CHR22_Formated1.gcount --bim ASD_XWASQC_chrX_sexFixed_Norm_Ref_Phased_EUR_Only${pop}_Info.bim -o chrX_sdMAF_ASD${pop}_XWAS1_CHR22_Formated -l LOG_chrX_sdMAF_ASD${pop}_XWAS1_CHR22_Formated ; done

