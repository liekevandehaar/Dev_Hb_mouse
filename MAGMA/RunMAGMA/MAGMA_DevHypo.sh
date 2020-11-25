#MAGMA for subsampled developmental Hb dataset

#!/bin/bash

#BIP
magma --gene-results ../../../../data/output/GWAS/BIP/G-BP_GWAS_Hou_2016.genes.raw \
    --gene-covar ../../../../data/input/Nov2020_GWAS_DevHyposamp_ANOVA.txt \
    --model direction=greater --out ../../../../data/output/GWAS/BIP/Nov2020_GWAS_BIP_DevHypo_ANOVA

#BMI
magma --gene-results ../../../../data/output/GWAS/BMI/G-BMI_GWAS_Locke_2015.genes.raw \
    --gene-covar ../../../../data/input/Nov2020_GWAS_DevHyposamp_ANOVA.txt \
    --model direction=greater --out ../../../../data/output/GWAS/BMI/Nov2020_GWAS_BMI_DevHypo_ANOVA

#Height
magma --gene-results ../../../../data/output/GWAS/Height/G-Height_GWAS_Wood_2015.genes.raw \
    --gene-covar ../../../../data/input/Nov2020_GWAS_DevHyposamp_ANOVA.txt \
    --model direction=greater --out ../../../../data/output/GWAS/Height/Nov2020_GWAS_Height_DevHypo_ANOVA

#MDD
magma --gene-results ../../../../data/output/GWAS/MDD2-18/G-MDD2_GWAS_2018.genes.raw \
    --gene-covar ../../../../data/input/Nov2020_GWAS_DevHyposamp_ANOVA.txt \
    --model direction=greater --out ../../../../data/output/GWAS/MDD2-18/Nov2020_GWAS_MDD2-18_DevHypo_ANOVA

#SCZ
magma --gene-results ../../../../data/output/GWAS/SCZ/G-SCZ_GWAS_Ripke_2014.genes.raw \
    --gene-covar ../../../../data/input/Nov2020_GWAS_DevHyposamp_ANOVA.txt \
    --model direction=greater --out ../../../../data/output/GWAS/SCZ/Nov2020_GWAS_SCZ_DevHypo_ANOVA
