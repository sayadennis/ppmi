import os
import numpy as np
import pandas as pd
from datetime import datetime

din = "/share/fsmresfiles/pd_project/wes/annovar_functionalfilter_new"
dout = "/share/fsmresfiles/pd_project/wes/variant_count_and_label"

mx_filename = "201007_ppmi_wes_variant_count_table.csv"

row_index = []
for fn in os.listdir(din):
    patidstring = fn[0:12] # PPMI_SI_<num>
    row_index.append(patidstring)

row_index.sort()

mx = pd.DataFrame(0, index=row_index, columns=[])

for i in range(len(row_index)):
    # First define ID of the row in mx
    patid = row_index[i]

    # Find the file name that corresponds to this row of mx
    for name in os.listdir(din):
        if row_index[i] in name:
            fn = name # full path or just filename? 
        else:
            continue
    
    # Read contents of filename into table, retrieve gene name and count variants
    # table = pd.read_csv(os.path.join(din, fn), sep="\t") ## CHANGE THIS --NO LONGER CSV. tab-delimited VCF 
    with open(os.path.join(din, fn), "r") as f:
        flines = f.readlines()
        temp_list = []
        for line in flines:
            if line.startswith("##"): # header
                continue
            elif line.startswith("#CHROM"): # colnames
                colnames = line.split("\t") 
            else:
                temp_list.append([line.split("\t")])
    
    for j in range(len(temp_list)): # read through each non-header line of VCF
        # get refGene's gene name
        row = temp_list[j][0]
        refgene = row[7].split(";")[np.where(["Gene.refGene" in x for x in row[7].split(";")])[0][0]]
        if "Gene.refGene" in refgene:
            refgene_val = refgene[13:]
            if "\\x3b" in refgene_val:
                refgene_split = refgene_val.split("\\x3b")
                refgene_concatenated = ""
                for l in range(len(refgene_split)-1):
                    refgene_concatenated = refgene_concatenated + refgene_split[l] + ";"
                refgene_concatenated = refgene_concatenated + refgene_split[-1]
                refgene_name = refgene_concatenated
            else:
                refgene_name = refgene_val
        else:
            print("Indexing failed to capture refGene gene name. Row: {}".format(row))
        # get variant count
        ct = 0
        gt_vals = row[9].split(":")[0].split("/")
        for k in gt_vals:
            if int(k) >= 1:
                ct += 1
            else:
                continue
        
        # write in VC table
        if refgene_name in list(mx.columns): # if there is already a column in mx named gene_name, put count value in mx.iloc[i, <column location>]
            mx.loc[patid][refgene_name] += ct
        else: # if not, create a new column named gene_names and put count calue in mx.iloc[i, <new column location>]
            mx[refgene_name] = 0 # add new column with default value zero 
            mx.loc[patid][refgene_name] = ct
    
    print("Done counting for patient {}... {}".format(patid, datetime.now().strftime("%m-%d-%Y, %H:%M:%S")))

mx.to_csv(path_or_buf=os.path.join(dout, mx_filename), header=True, index=True) 
