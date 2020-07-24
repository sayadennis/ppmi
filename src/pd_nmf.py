import numpy as np
from sklearn.decomposition import NMF

def optimize_k(k_list, X_train, y_train, X_test, y_test, hypertune_func, param_grid, multi_class): # X and y should be numpy arrays. Rows of X are patients.
    
    train_size = X_train.shape[0] # rows are patients for X
    X = np.concatenate((X_train, X_test))
    y = np.concatenate((y_train, y_test))
    
    performance_dict = {}
    
    for k in k_list:
        A = X.T
        nmf = NMF(n_components=k, init="nndsvd", random_state=24)
        W = nmf.fit_transform(A)
        F = np.dot(A.T, W)
        F_train = F[:train_size, :]
        F_test = F[train_size:, :]
        
        opt_params, opt_mean_score = hypertune_func(F_train, y_train, param_grid, multi=multi_class)
        
        performance_dict[k] = [opt_mean_score, opt_params]
    return F_test, performance_dict


def get_F(k, X_train, y_train, X_test, y_test):
    train_size = X_train.shape[0] # rows are patients for X
    X = np.concatenate((X_train, X_test))
    y = np.concatenate((y_train, y_test))

    A = X.T
    nmf = NMF(n_components=k, init="nndsvd", random_state=24)
    W = nmf.fit_transform(A)
    F = np.dot(A.T, W)
    F_train = F[:train_size, :]
    F_test = F[train_size:, :]

    return F_train, F_test, W


def coef_genes(W, genes, thres=0.01): # length of genes and rows of W should match
    if W.shape[0] != len(genes):
        return -1
    else:
        group_list = []
        for i in range(W.shape[1]): # iterate through columns of W i.e. weight vectors of each factor groups 
            coef_genes = genes[np.where(W[:,i] > thres)[0]]
            genes_final = []
            for j in range(len(coef_genes)):
                split_list = coef_genes[j].split(";") # element in list can contain multiple gene names (overlapping)
                for gene in split_list:
                    if gene in genes_final:
                        continue
                    else:
                        genes_final.append(gene)
            group_list.append([genes_final])
    return genes_final

