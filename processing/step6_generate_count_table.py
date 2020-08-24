# The purpose of this script is to parse through the PPMI_SI_<num>_exac.txt files
# Find the genotype count (0/1/2) and write a line with variant counts in the appropriate gene column into the final matrix.

import os
import pandas as pd
from datetime import datetime

input_path = "/share/fsmresfiles/pd_project/wes/filtered_exac"
output_path = "/share/fsmresfiles/pd_project/wes/variant_count_and_label"

mx_filename = "200824_ppmi_wes_variant_count_table.csv"

row_index = []
for filename in os.listdir(input_path):
    string = filename[0:12]
    row_index.append(string)

row_index.sort()

mx = pd.DataFrame(0, index=row_index, columns=[])

for i in range(len(row_index)):
    # First define ID of the row in mx
    patient_id = row_index[i]

    # Find the file name that corresponds to this row of mx
    for name in os.listdir(input_path):
        if row_index[i] in name:
            filename = name
        else:
            continue
    
    # Read contents of filename into table, retrieve gene name and count variants
    table = pd.read_csv(os.path.join(input_path, filename), sep="\t")
    for j in range(table.shape[0]): # loop through lines of table
        gene_name = table.iloc[j,6] # refGene
        ct = 0
        gt = table.iloc[j,77].split(":")[0] # gt = "0/0", "1/0", "2/2" etc.
        gt_values = [int(gt.split("/")[0]), int(gt.split("/")[1])]
        for k in gt_values:
            if k >= 1:
                ct += 1
            else:
                continue
        
        # if there is already a column in mx named gene_name, put count value in mx.iloc[i, <column location>]
        # if not, create a new column named gene_names and put count calue in mx.iloc[i, <new column location>]
        if gene_name in list(mx.columns):
            mx.iloc[i, mx.columns==gene_name] += ct
        else:
            mx[gene_name] = 0 # add new column with default value zero 
            mx.iloc[i, mx.columns==gene_name] = ct
    
    print("Done counting for patient {}... {}".format(patient_id, datetime.now().strftime("%m-%d-%Y, %H:%M:%S")))

mx.to_csv(path_or_buf=os.path.join(output_path, mx_filename), sep="\t")
