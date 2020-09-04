import os
import glob
import pandas as pd
from datetime import datetime

input_path = "/share/fsmresfiles/pd_project/wes/annovar_vcf"
output_path = "/share/fsmresfiles/pd_project/wes/variant_count_and_label"

fn_allmx = "200825_ppmi_wes_variant_count_table_include_nonLGD.csv"
fn_snpmx = "200825_snp_wes_vc_include_nonLGD.csv"
fn_indelmx = "200825_indel_wes_vc_include_nonLGD.csv"


row_index = []
for pn in glob.glob(os.path.join(input_path, "*_anv.hg19_multianno.txt")): # pathname
    row_index.append(pn[46:58])

row_index.sort()

allmx = pd.DataFrame(0, index=row_index, columns=[])
snpmx = pd.DataFrame(0, index=row_index, columns=[])
indelmx = pd.DataFrame(0, index=row_index, columns=[])


for i in range(len(row_index)):
    # First define ID of the row in mx
    patient_id = row_index[i]

    # Find the file name that corresponds to this row of mx
    for pn in glob.glob(os.path.join(input_path, "*_anv.hg19_multianno.txt")): # pathname
        if row_index[i] in pn:
            filepath = pn
        else:
            continue
    
    # Read contents of filename into table, retrieve gene name and count variants
    table = pd.read_csv(filepath, sep="\t", skiprows=[0], header=None)
    colnames = [
        "CHROM", "POS", "ID", "REF", "ALT", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene", 
        "Func.knownGene", "Gene.knownGene", "GeneDetail.knownGene", "ExonicFunc.knownGene", "AAChange.knownGene", 
        "Func.ensGene", "Gene.ensGene", "GeneDetail.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene", 
        "snp138", "SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred", 
        "LRT_score", "LRT_pred", "MutationTaster_score", "MutationTaster_pred", "MutationAssessor_score", "MutationAssessor_pred", 
        "FATHMM_score", "FATHMM_pred", "PROVEAN_score", "PROVEAN_pred", "VEST3_score", "CADD_raw", "CADD_phred", "DANN_score", 
        "fathmm-MKL_coding_score", "fathmm-MKL_coding_pred", "MetaSVM_score", "MetaSVM_pred", "MetaLR_score", "MetaLR_pred", 
        "integrated_fitCons_score", "integrated_confidence_value", "GERP++_RS", "phyloP7way_vertebrate", "phyloP20way_mammalian", 
        "phastCons7way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds", "Interpro_domain", "1000g2015aug_all", 
        "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", 
        "GT", "QUAL", "DP", "chr", "pos", "name", "ref", "alt", "qual", "FILTER", "INFO", "FORMAT", "FORMAT_INFO"
    ]
    table.columns = colnames
    for j in range(table.shape[0]): # loop through lines of table
        gene_name = table.iloc[j,6] # refGene
        all_ct = 0
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
        
        for (mx, ct) in [(allmx, all_ct), (snpmx, snp_ct), (indelmx, indel_ct)]:
        # if there is already a column in mx named gene_name, put count value in mx.iloc[i, <column location>]
        # if not, create a new column named gene_names and put count calue in mx.iloc[i, <new column location>]
            if gene_name in list(mx.columns):
                mx.iloc[i, mx.columns==gene_name] += ct # add current count to the existing count
            else:
                mx[gene_name] = 0 # add new column with default value zero. Is this correct?
                mx.iloc[i, mx.columns==gene_name] = ct
    
    print("Done counting for patient {}... {}".format(patient_id, datetime.now().strftime("%m-%d-%Y, %H:%M:%S")))

allmx.to_csv(os.path.join(output_path, fn_allmx), header=True, index=True)
snpmx.to_csv(os.path.join(output_path, fn_snpmx), header=True, index=True)
indelmx.to_csv(os.path.join(output_path, fn_indelmx), header=True, index=True)
