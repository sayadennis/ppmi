# The purpose of this script is to build the label numpy array for ppmi wes variant count matrix

import os
import pandas as pd

input_path = "../data/phenotype"
wes_path = "../data/wes/filtered_exac"
output_path = "../data/wes"

input_filename = "primary_diagnosis.csv"

input_table = pd.read_csv(os.path.join(input_path, input_filename), sep=",")

## Generate row index

row_index = []
for filename in os.listdir(wes_path):
    string = filename[0:12]
    row_index.append(string)

row_index.sort()




## Later timepoint

output_filename = "191113_primdiag_latertimepoint.csv"

output_table = pd.DataFrame(index=row_index, columns=["Label"])
output_table = output_table[output_table.index != "PPMI_SI_3329"]

for patient_id in output_table.index:
    label = input_table.loc[input_table["PATNO"]==int(patient_id[-4:]), "PRIMDIAG"]
    output_table.loc[patient_id]["Label"] = int(label.values[-1])

# output_table.to_csv(path_or_buf=os.path.join(output_path, output_filename), sep="\t")

## PD at some point: 1 is PD

output_filename = "191113_primdiag_pdanypoint.csv"

output_table = pd.DataFrame(index=row_index, columns=["Label"])
output_table = output_table[output_table.index != "PPMI_SI_3329"]

for patient_id in output_table.index:
    label = input_table.loc[input_table["PATNO"]==int(patient_id[-4:]), "PRIMDIAG"]
    if any(label.values==1.):
        output_table.loc[patient_id]["Label"] = 1
    else:
        output_table.loc[patient_id]["Label"] = int(label.values[-1])


# output_table.to_csv(path_or_buf=os.path.join(output_path, output_filename), sep="\t")

## Building binary of HC (0) vs anything else (1) --410 "1"s.

output_filename = "191113_primdiag_binary_hc_or_not.csv"

for patient_id in output_table.index:
    if output_table.loc[patient_id]["Label"] == 17:
        output_table.loc[patient_id]["Label"] = 0
    else:
        output_table.loc[patient_id]["Label"] = 1

# output_table.to_csv(path_or_buf=os.path.join(output_path, output_filename), sep="\t")

## Building a binary of PD (1) vs anything else (0) --396 "1"s.

output_filename = "191113_primdiag_binary_pd_or_not.csv"

for patient_id in output_table.index:
    if output_table.loc[patient_id]["Label"] == 1:
        output_table.loc[patient_id]["Label"] = 1
    else:
        output_table.loc[patient_id]["Label"] = 0

# output_table.to_csv(path_or_buf=os.path.join(output_path, output_filename), sep="\t")
