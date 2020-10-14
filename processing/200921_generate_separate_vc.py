import os
import numpy as np
import pandas as pd
from datetime import datetime

din = "/share/fsmresfiles/pd_project/wes/annovar_qualityfilter_VariantAnnotation"
dout = "/share/fsmresfiles/pd_project/wes/variant_count_and_label"

exonicfunc_list = [
    "nonsynonymous_SNV", 
    "synonymous_SNV", 
    "nonframeshift_insertion", 
    "unknown", 
    "stopgain", 
    "nonframeshift_deletion", 
    "frameshift_insertion", 
    "frameshift_deletion", 
    "stoploss",
    "missing"
]

func_list = [
    "upstream",
    "ncRNA_exonic",
    "downstream",
    "intergenic",
    "ncRNA_intronic",
    "intronic",
    "exonic",
    "UTR3",
    "UTR5",
    "upstream_downstream",
    "splicing",
    "ncRNA_splicing",
    "UTR5_UTR3",
    "exonic_splicing",
    "ncRNA_exonic_splicing",
    "missing"
]

fn_dict = {}
for i in range(len(exonicfunc_list)):
    fn_dict[exonicfunc_list[i]] = "200921_vc_{}.csv".format(exonicfunc_list[i])


row_index = []
for fn in os.listdir(din):
    patidstring = fn[0:12] # PPMI_SI_<num>
    row_index.append(patidstring)

row_index.sort()

mx_dict = {}
for i in range(len(exonicfunc_list)):
    mx_dict[exonicfunc_list[i]] = pd.DataFrame(None, index=row_index, columns=[])


## Now, we loop through the files based on row indices and fill in the matrices 

for i in range(len(row_index)):
    # First define ID of the row in mx
    patid = row_index[i]

    # Find the file name that corresponds to this row of mx
    for name in os.listdir(din):
        if row_index[i] in name:
            fn = name # just filename, not full path
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
        
        # get exonic function
        exonicfunc = row[7].split(";")[np.where(["ExonicFunc.refGene" in x for x in row[7].split(";")])[0][0]].split("=")[1]
        if exonicfunc == ".":
            exonicfunc = "missing"

        # get variant count
        ct = 0
        gt_vals = row[9].split(":")[0].split("/")
        for k in gt_vals:
            if int(k) >= 1:
                ct += 1
            else:
                continue
        
        # all matrices in mx_dict will have synchronized columns, for easier handling later on.

        # write in VC table
        if refgene_name in list(mx_dict["nonsynonymous_SNV"].columns): # if there is already a column in mx named gene_name, put count value in mx.iloc[i, <column location>]
            try:
                mx_dict[exonicfunc].loc[patid][refgene_name] += ct
            except:
                print("Error: Could not write to matrices –– Gene: {} // ExonicFunc: {} // Count: {}".format(refgene_name, exonicfunc, ct))
        else: # if not, create a new column named gene_names and put count calue in mx.iloc[i, <new column location>]
            for key in mx_dict:
                mx_dict[key][refgene_name] = 0 # add new column with default value zero to ALL matrices
            mx_dict[exonicfunc].loc[patid][refgene_name] = ct # add value to the matrix with corresponding vartype 
    
    print("Done counting for patient {}... {}".format(patid, datetime.now().strftime("%m-%d-%Y, %H:%M:%S")))

for key in mx_dict:
    mx_dict[key].to_csv(path_or_buf=os.path.join(dout, fn_dict[key]), header=True, index=True)
