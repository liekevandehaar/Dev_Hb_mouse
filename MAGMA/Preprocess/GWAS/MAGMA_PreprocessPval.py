##### Format P-value from scientific format to decimal format #####
# author: Juliska E Boer
# date: 04 Nov 2020

#!/bin/python

#the p-value is in column 9, and is converted to decimal format with 10 decimals
with open("Height_GWAS_format.txt", "r") as f, open("../../data/input/Height_GWAS_final.txt", "w") as writef:
    for lines in f:
        line = lines.split("\t")
        new_line = line
        if "e" in line[8]:
            new_line[8] = format(float(line[8]), ".10f")
        writef.writelines('\t'.join(new_line))