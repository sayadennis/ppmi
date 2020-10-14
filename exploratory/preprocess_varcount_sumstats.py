import os
import numpy as np
import pandas as pd

pre_fn = "200922_varcounts_prefilter.out"
qual_fn = "200922_varcounts_qualfilter.out"
func_fn = "200922_varcounts_funcfilter.out"

pre_mx = pd.DataFrame(None, columns=["SNPs", "Indels", "Total"])
qual_mx = pd.DataFrame(None, columns=["SNPs", "Indels", "Total"])
func_mx = pd.DataFrame(None, columns=["SNPs", "Indels", "Total"])

with open(pre_fn) as f:
    lines = f.readlines()
    # There is 7 lines for each patient 
    indices = np.where([line.startswith("##") for line in lines])[0]
    for i in range(len(indices)): # indices of where lines start with "##"
        try:
            patid = lines[indices[i]][2:-1] # since last character is \n
            snps = int(lines[indices[i]+2].split()[-1])
            indels = int(lines[indices[i]+3].split()[-1])
            total = int(lines[indices[i]+6].split()[-1])
            pre_mx = pre_mx.append(pd.DataFrame(np.array([[snps, indels, total]]), index=[patid], columns=["SNPs", "Indels", "Total"]))
        except: # if it gives an error, ignore it (3332 and 3504 or something)
            continue


with open(qual_fn) as f:
    lines = f.readlines()
    # There is 7 lines for each patient 
    indices = np.where([line.startswith("##") for line in lines])[0]
    for i in range(len(indices)): # indices of where lines start with "##"
        try:
            patid = lines[indices[i]][2:-1] # since last character is \n
            snps = int(lines[indices[i]+2].split()[-1])
            indels = int(lines[indices[i]+3].split()[-1])
            total = int(lines[indices[i]+6].split()[-1])
            qual_mx = qual_mx.append(pd.DataFrame(np.array([[snps, indels, total]]), index=[patid], columns=["SNPs", "Indels", "Total"]))
        except: # if it gives an error, ignore it (3332 and 3504 or something)
            continue


with open(func_fn) as f:
    lines = f.readlines()
    # There is 7 lines for each patient 
    indices = np.where([line.startswith("##") for line in lines])[0]
    for i in range(len(indices)): # indices of where lines start with "##"
        try:
            patid = lines[indices[i]][2:-1] # since last character is \n
            snps = int(lines[indices[i]+2].split()[-1])
            indels = int(lines[indices[i]+3].split()[-1])
            total = int(lines[indices[i]+6].split()[-1])
            func_mx = func_mx.append(pd.DataFrame(np.array([[snps, indels, total]]), index=[patid], columns=["SNPs", "Indels", "Total"]))
        except: # if it gives an error, ignore it (3332 and 3504 or something)
            continue

pre_mx.to_csv("/home/srd6051/pd_project/prefilter_varcounts.csv", header=True, index=True)
qual_mx.to_csv("/home/srd6051/pd_project/qualfilter_varcounts.csv", header=True, index=True)
func_mx.to_csv("/home/srd6051/pd_project/funcfilter_varcounts.csv", header=True, index=True)
