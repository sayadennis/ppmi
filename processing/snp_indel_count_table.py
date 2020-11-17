import os
from datetime import datetime
import pandas as pd

input_path = "/data/filtered_exac"
output_path = "/data/variant_count_and_label"

fn_snpmx = "snp_wes_vc.csv"
fn_indelmx = "indel_wes_vc.csv"

# generate row index as list
row_index = []
for filename in os.listdir(input_path):
    string = filename[0:12]
    row_index.append(string)

row_index.sort()

snpmx = pd.DataFrame(0, index=row_index, columns=[])
indelmx = pd.DataFrame(0, index=row_index, columns=[])

for i in range(len(row_index)):
    # First define patient ID
    patient_id = row_index[i]

    # Find the file name that corresponds to this row of mx
    for name in os.listdir(input_path):
        if row_index[i] in name:
            fn = name
        else:
            continue
    
    # Read contents of fn into table, retrieve gene name and count variants
    table = pd.read_csv(os.path.join(input_path, fn), sep="\t", header=0, index_col=0)
    for j in range(table.shape[0]): # loop through lines of table
        gene_name = table.iloc[j,6] # refGene gene name
        snp_ct = 0
        indel_ct = 0
        gt = table.iloc[j,77].split(":")[0] # gt = "0/0", "1/0", "2/2" etc.
        gt_values = [int(gt.split("/")[0]), int(gt.split("/")[1])]
        if len(table.iloc[j,3]) != len(table.iloc[j,4]): ## FIX THIS PART
            for k in gt_values:
                if k >= 1:
                    indel_ct += 1
                else:
                    continue
        elif table.iloc[j,4] == ".": ## ALSO REVIEW THESE CRITERIA TO ENSURE I'M COUNTING CORRECTLY
            for k in gt_values:
                if k >= 1:
                    indel_ct += 1
                else:
                    continue
        else: ## THIS TOO
            for k in gt_values:
                if k >= 1:
                    snp_ct += 1
                else:
                    continue
        
        for (mx, ct) in [(snpmx, snp_ct), (indelmx, indel_ct)]:
            # if there is already a column in mx named gene_name, put count value in mx.iloc[i, <column location>]
            # if not, create a new column named gene_names and put count calue in mx.iloc[i, <new column location>]
            if gene_name in list(mx.columns):
                mx.iloc[i, mx.columns==gene_name] += ct # add current count to the existing count
            else:
                mx[gene_name] = 0 # add new column with default value zero. Is this correct?
                mx.iloc[i, mx.columns==gene_name] = ct
        
    print("Done counting for patient {}... {}".format(patient_id, datetime.now().strftime("%m-%d-%Y, %H:%M:%S")))


snpmx.to_csv(path_or_buf=os.path.join(output_path, fn_snpmx), header=True, index=True)
indelmx.to_csv(path_or_buf=os.path.join(output_path, fn_indelmx), header=True, index=True)
