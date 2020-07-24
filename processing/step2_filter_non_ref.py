# The purpose of this script is to filter the VCF files by: 
    # loop through the raw VCFs
    # remove lines where only "<NON_REF>" is the alternative allele 
    # and also remove lines that reference "<NON_REF>" in the GT column

import os

rawvcf_dir = "../data/wes/raw_vcf"
filtered_dir = "../data/wes/nonref_filtered"

for oldfilename in os.listdir(rawvcf_dir):
    if oldfilename.endswith(".vcf"):
        with open(os.path.join(rawvcf_dir, oldfilename), "r") as old_vcf:
            lines = old_vcf.readlines()
            newfilename = oldfilename.split(".")[0] + "_filter.vcf"
        with open(os.path.join(filtered_dir, newfilename), "w") as new_vcf:
            for line in lines:
                if "\t<NON_REF>\t" in line: # cases where only alternative allele is <NON_REF>
                    continue
                elif ",<NON_REF>" in line: # cases where alternative allele includes <NON_REF> and some other allele
                    cols = line.split("\t")
                    alt = cols[4] # alternative allele is the fourth column (double check!!)
                    alt_lst = alt.split(",") # list of alternative alleles
                    ind = cols[9] # individual's information
                    ind_gt = ind.split(":") # the GT of the individual
                    gt_num = ind_gt[0].split("/") # the first part of the GT, which references the alt allele 
                    # now alt_lst contains the alternative alleles - list of size 1 or larger
                    # and gt_num contains a list of two numbers (in string form) of which alternative allele that line is referencing
                    if max(int(gt_num[0]), int(gt_num[1])) > len(alt_lst): # cases where GT references the alt allele
                        continue
                    else:
                        sub = line.split(",<NON_REF>")
                        new_vcf.write(sub[0] + sub[1])
                else:
                    new_vcf.write(line)
    else: # ignore files that end with .idx
        continue
