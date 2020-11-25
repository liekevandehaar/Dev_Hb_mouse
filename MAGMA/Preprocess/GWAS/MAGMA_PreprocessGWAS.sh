##### Preprocessing GWAS data for MAGMA analysis #####
# author: Juliska E Boer
# date: 04 Nov 2020

#!/bin/bash

#use this file for formatting
awk -F "\t" '{print $2"\t"$1":"$4}' ../../data/input/g1000_eur.bim > g1000_formatting.bim

#format BMI data
#misses the CHR and BP column

#retrieve CHR and BP through SNP ID from formatting file
merged=$(join -1 1 -2 2 -o 1.1,2.1,2.4,1.2,1.3,1.4,1.5,1.6,1.7,1.8 <(sort -k1 ../../data/input/BMI_GWAS_Locke_2015.uniq) <(sort -k2 ../../data/input/g1000_eur.bim) | tr " " "\t")
( echo -e "SNP\tCHR\tBP\tA1\tA2\tFreq1.Hapmap\tb\tse\tP\tN"; echo -e "${merged}\n" ) > ../../data/input/BMI_GWAS_final.txt

#format BIP data
#misses the SNP ID column, has MarkerName instead

#change order of columns
awk -F "\t" '{split($3,a,":"); print(a[1]":"a[2]"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17)}' ../../data/input/BP_GWAS_Hou_2016.txt > BP_GWAS_order.txt
#retrieve SNP ID through MarkerName and formatting file
(join -1 1 -2 2 -o 2.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17 <(sort -k1 BP_GWAS_order.txt) <(sort -k2 g1000_formatting.bim)) | tr " " "\t" > BP_GWAS_format.txt
#add header
header="SNP\tCHR\tBP\tAllele1\tAllele2\tFreq1\tFreqSE\tWeight\tZscore\tP\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\tEffect\tStdErr"
echo -e ${header} > outfile ; cat BP_GWAS_format.txt >> outfile
mv outfile ../../data/input/BP_GWAS_final.txt
rm BP_GWAS_order.txt
rm BP_GWAS_format.txt

#format Height data
#misses the CHR and BP column
#p-value incorrect format (scientific)

#retrieve CHR and BP through SNP ID from formatting file
(join -1 1 -2 2 -o 1.1,2.1,2.4,1.2,1.3,1.4,1.5,1.6,1.7,1.8 <(sort -k1 ../../data/input/Height_GWAS_Wood_2015.txt) <(sort -k2 ../../data/input/g1000_eur.bim) | tr " " "\t") > Height_GWAS_snps.txt
#add header
header="SNP\tCHR\tBP\tAllele1\tAllele2\tFreq.Hapmap\tB\tSE\tP\tN"
echo -e ${header} > outfile ; cat Height_GWAS_snps.txt >> outfile
mv outfile Height_GWAS_format.txt
rm Height_GWAS_snps.txt
#p-value: scientific format to decimal format
python3 FormatPvals.py
rm Height_GWAS_format.txt

#format Major Depressive Disorder
#the order of the columns is not correct

#reorder columns
cat ../../data/input/MDD2018_ex23andMe.txt | sed '/ /d' | awk -F "\t" '{print ($2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17)}' > ../../data/input/MDD2_GWAS_final.txt

#format Schizophrenia data
#order of columns is not correct

awk -F "\t" '{print ($2"\t"$1"\t"$5"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10)}' ../../data/input/SCZ_GWAS_Ripke_2014.scz2snpres | awk  'BEGIN {OFS="\t"} {$2 = substr($2, 4); print}' > SCZ_GWAS_format.txt
#add header
header="SNP\tCHR\tBP\tAllele1\tAllele2\tInfo\tOR\tSE\tP\tNGT"
echo -e ${header} > outfile ; cat SCZ_GWAS_format.txt | tail -n +2 >> outfile
mv outfile ../../data/input/SCZ_GWAS_final.txt
rm SCZ_GWAS_format.txt

rm g1000_formatting.bim


