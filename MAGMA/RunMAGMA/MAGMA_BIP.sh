##### MAGMA analysis for Height GWAS #####
# author: Juliska E Boer
# date: 04 Nov 2020

#!/bin/bash

#perform gene annotation
magma --annotate window=1 --snp-loc ../../data/input/BP_GWAS_final.txt \
    --gene-loc ../../data/input/ENSGv92.coding.genes.txt \
    --out ../../data/output/GWAS/BIP/A-BP_GWAS_Hou_2016

#perform gene analysis
magma --bfile ../../data/input/g1000_eur/g1000_eur \
    --pval ../../data/input/BP_GWAS_final.txt use=1,10 N=34950 \
    --gene-annot ../../data/output/GWAS/BIP/A-BP_GWAS_Hou_2016.genes.annot \
    --out ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016

#gene property analysis: protocol 1 by Nathan Skene
#Developmental Hb - protocol 1
magma --gene-results ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../data/input/Aug2020_GWAS_embryo_ANOVA.txt \
    --model direction=greater --out ../../data/output/GWAS/BIP/Aug2020_GWAS_BIP_embryo_ANOVA
  
#Hashikawa, et al (2020) - protocol 1
magma --gene-results ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../data/input/Aug2020_GWAS_stuber_ANOVA.txt \
    --model direction=greater --out ../../data/output/GWAS/BIP/Aug2020_GWAS_BIP_stuber_ANOVA
  
#Wallace, et al (2020) - protocol 1
magma --gene-results ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../data/input/Aug2020_GWAS_wallace_ANOVA.txt \
    --model direction=greater --out ../../data/output/GWAS/BIP/Aug2020_GWAS_BIP_wallace_ANOVA
	
#Tabula Muris, et al (2018) - protocol 1
magma --gene-results ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../data/input/Aug2020_GWAS_tabulam_ANOVA.txt \
    --model direction=greater --out ../../data/output/GWAS/BIP/Aug2020_GWAS_BIP_TabulaMuris_ANOVA

#Merged mouse Hb - protocol 1
magma --gene-results ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../data/input/Oct2020_GWAS_merged_ANOVA.txt \
    --model direction=greater --out ../../data/output/GWAS/BIP/Oct2020_GWAS_BIP_merged_ANOVA
	
#gene property analysis: protocol 2 by Kyoko Watanabe
#Developmental Hb - protocol 2
magma --gene-results ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../data/input/Aug2020_GWAS_embryo_AVG.txt \
    --model direction=greater condition-hide=AVERAGE --out ../../data/output/GWAS/BIP/Aug2020_GWAS_BIP_embryo_AVG
  
#Hashikawa, et al (2020) - protocol 2
magma --gene-results ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../data/input/Aug2020_GWAS_stuber_AVG.txt \
    --model direction=greater condition-hide=AVERAGE --out ../../data/output/GWAS/BIP/Aug2020_GWAS_BIP_stuber_AVG
  
#Wallace, et al (2020) - protocol 2
magma --gene-results ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../data/input/Aug2020_GWAS_wallace_AVG.txt \
    --model direction=greater condition-hide=AVERAGE --out ../../data/output/GWAS/BIP/Aug2020_GWAS_BIP_wallace_AVG

#Tabula Muris, et al (2018) - protocol 2
magma --gene-results ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../data/input/Aug2020_GWAS_tabulam_AVG.txt \
    --model direction=greater condition-hide=AVERAGE --out ../../data/output/GWAS/BIP/Aug2020_GWAS_BIP_TabulaMuris_AVG

#Merged mouse Hb - protocol 2
magma --gene-results ../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../data/input/Oct2020_GWAS_merged_AVG.txt \
    --model direction=greater condition-hide=AVERAGE --out ../../data/output/GWAS/BIP/Oct2020_GWAS_BIP_merged_AVG