import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

def odds_ratio(X, y, n_classes=3):
    odds_array = np.zeros((X.shape[1], n_classes)) # there are 3 subtypes
    for i in range(X.shape[1]):
        outcomes_present = y[X[:,i] > 0]
        outcomes_absent = y[X[:,i] == 0]
        for j in range(n_classes):
            num = (outcomes_present == np.unique(y)[j]).sum() / len(outcomes_present)
            denom = (outcomes_absent == np.unique(y)[j]).sum() / len(outcomes_absent)
            odds_array[i,j] = num/denom
    return odds_array

# function that takes input and target, returns a list of gene names that had a significant p-value
def or_sig_genes(X_table, y_table, thr=0.05, n_classes=3): # X and y are pandas dataframes
    X = X_table.to_numpy()
    y = y_table.to_numpy().reshape((-1,))
    table_pvals = np.empty((X.shape[1], n_classes)) # n_classes = different outcomes
    table_pvals[:] = np.NaN
    for i in range(X.shape[1]): # iterate through each gene
        outcomes_present = y[X[:,i] > 0]
        outcomes_absent = y[X[:,i] == 0]
        for j in range(n_classes):
            present_positive = (outcomes_present == j+1).sum()
            present_negative = (outcomes_present != j+1).sum()
            absent_positive = (outcomes_absent == j+1).sum()
            absent_negative = (outcomes_absent != j+1).sum()
            contingency_table = np.array([[present_positive, present_negative], [absent_positive, absent_negative]])
            _, pval = fisher_exact(contingency_table, alternative="two-sided")
            table_pvals[i,j] = pval
    genes = X_table.columns[np.where(table_pvals < thr)[0]] # this list contains duplicates
    genes_final = []
    for i in range(len(genes)):
        split_list = genes[i].split(";") # element in list can contain multiple gene names (overlapping)
        for gene in split_list:
            if gene in genes_final:
                continue
            else:
                genes_final.append(gene)
    return table_pvals, genes_final
