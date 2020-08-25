import numpy as np
import pandas as pd
from sklearn import metrics

from sklearn.linear_model import LogisticRegressionCV
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import OneHotEncoder

from sklearn.decomposition import NMF

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif
from sklearn.feature_selection import mutual_info_classif
from sklearn.feature_selection import chi2
from sklearn.feature_selection import f_regression


class HyperparameterTuning:
    """
    This class facilitates custom hyperparameter tuning.

    Usage:
        ht = HyperparameterTuning(X, y)
        opt_params, opt_mean_score = ht.tune_lrm(param_list)
        opt_params, opt_mean_score = ht.tune_svm(param_list)
    """

    def __init__(self, multi=True):
        self.multi = multi
        self.lrm_Cs = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1., 1e+1, 1e+2, 1e+3]
        self.svm_params = {
            "kernel" : ["linear", "rbf"], 
            "C" : [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e+0, 1e+1, 1e+2, 1e+3], 
            "gamma" : [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1., 1e+1, 1e+2, 1e+3]
            }
    
    def tune_lrm(self, X, y, param_list=None):
        # set default parameter grid
        if param_list == None:
            param_list = self.lrm_Cs
        
        # specify scoring and refit criteria respectively for multiclass and binary
        if self.multi:
            multi_class = "multinomial"
            scoring = "roc_auc_ovr_weighted"
        else:
            multi_class = "auto"
            scoring = "roc_auc"
        
        # cross-validate best parameter 
        lrCV = LogisticRegressionCV(Cs=param_list, class_weight="balanced", n_jobs=8, max_iter=3000,
            multi_class=multi_class, scoring=scoring, cv=StratifiedKFold(n_splits=5, random_state=1, shuffle=True), refit=True
        )
        lrCV.fit(X, y)
        opt_C = lrCV.C_[0]
        opt_mean_score = np.mean([lrCV.scores_[x][:,np.where(param_list == opt_C)] for x in range(1,4)]) # for classes 1-3 

        return opt_C, opt_mean_score
    
    def tune_svm(self, X, y, param_list=None):
        # set default parameter grid
        if param_list == None:
            param_list = self.svm_params

        # specify scoring and refit criteria respectively for multiclass and binary
        if self.multi:
            scoring = "roc_auc_ovr_weighted"
            refit = True
        else:
            scoring = "roc_auc"
            refit = True

        # gridsearch and cross-validate best set of parameters 
        gsCV = GridSearchCV(
            SVC(class_weight="balanced", max_iter=3000, probability=True, random_state=240, decision_function_shape="ovr"), 
            param_grid=param_list, n_jobs=8, scoring=scoring, refit=refit,
            cv=StratifiedKFold(n_splits=5, random_state=1, shuffle=True)
        )
        gsCV.fit(X, y)
        opt_params = gsCV.best_params_
        opt_kernel = gsCV.best_params_["kernel"]
        opt_C = gsCV.best_params_["C"]
        opt_gamma = gsCV.best_params_["gamma"]
        opt_mean_score = np.mean(
            gsCV.cv_results_["mean_test_score"][
                (gsCV.cv_results_["param_kernel"] == opt_kernel) & 
                (gsCV.cv_results_["param_C"] == opt_C) & 
                (gsCV.cv_results_["param_gamma"] == opt_gamma)
            ]
        )

        return opt_params, opt_mean_score



class PD_NMF:
    """
    This class handles the processes surrounding NMF of PD variant data.

    Usage:
        nmf = PD_NMF(X_train, y_train)
        nmf.tune_k(k_list, model)
        print("Best performance {} with k={} with parameters {}".format(nmf.opt_mean_score, nmf.opt_k, nmf.opt_params))
        F_test = nmf.transform(X_test)
        pred = model(F_test) etc. etc. 
    """
    def __init__(self, X, y):
        self.X = X
        self.y = y
    
    def tune_k(self, k_list, hypertune_func, multi):
        performance_dict = {}
        for k in k_list:
            A = self.X.T
            nmf = NMF(n_components=k, init="nndsvd", random_state=24)
            W = nmf.fit_transform(A)
            F_train = np.dot(A.T, W)
            
            opt_params, opt_mean_score = hypertune_func(F_train, self.y) # multi=multi
            performance_dict[k] = [opt_mean_score, opt_params, W, F_train]

        self.opt_k = list(performance_dict.keys())[np.argmax([list(performance_dict.values())[i][0] for i in range(len(performance_dict))])]
        self.opt_mean_score = performance_dict[self.opt_k][0]
        self.opt_params = performance_dict[self.opt_k][1]
        self.opt_W = performance_dict[self.opt_k][2]
        self.opt_F = performance_dict[self.opt_k][3]
    
    def transform(self, X_test):
        return np.dot(X_test, self.opt_W)



class TuneSelectKBest(HyperparameterTuning):
    """
    This class handles the tuning of the parameter K as well as selection criteria for SelectKBest function.

    Example usage:
        ktuner = TuneSelectKBest(X_train, X_test, y_train, y_test)
        opt_k, opt_criteria, opt_C, roc_auc = ktuner.tune_selectkbest(k_list)
        gene_set = ktuner.extract_genes(vc.columns, n_genes=opt_k)
    """

    def __init__(self, X_train, X_test, y_train, y_test):
        super().__init__(multi=True) # inherit all methods from HyperparameterTuning class 
        self.X_train = X_train
        self.X_test = X_test
        self.y_train = y_train
        self.y_test = y_test
        self.criteria_dict = {
            "f_classif" : f_classif, 
            "mutual_info_classif" : mutual_info_classif, 
            "chi2" : chi2, 
        }
        self.skb = None # will later be defined when tune_selectkbest is called 
    
    def tune_selectkbest(self, k_list, multi=True):
        pf = np.zeros((len(k_list), len(self.criteria_dict.keys())))
        for k in k_list:
            for key in self.criteria_dict.keys():
                skb = SelectKBest(self.criteria_dict[key], k=k)
                skb.fit(self.X_train, self.y_train.astype(np.float).ravel())
                select_X_train = skb.transform(self.X_train)
                select_X_test = skb.transform(self.X_test)

                encoder = OneHotEncoder()
                y_onehot = encoder.fit_transform(self.y_test.reshape((-1,1))).toarray()

                # tune C for LRM
                ht = HyperparameterTuning()
                opt_C, _ = ht.tune_lrm(select_X_train, self.y_train.astype(np.float).ravel())

                # train LRM with optimal C
                lrm = LogisticRegression(C=opt_C, class_weight="balanced", max_iter=3000, n_jobs=6)
                lrm.fit(select_X_train, self.y_train.astype(np.float).ravel())
                
                # record performance
                y_decf = lrm.decision_function(select_X_test)

                # compute micro-averaged ROC-AUC
                fpr, tpr, _ = metrics.roc_curve(y_onehot.ravel(), y_decf.ravel())
                roc_auc = metrics.auc(fpr, tpr)

                pf[k_list == k, self.criteria_dict.keys() == key] = roc_auc
        
        opt_k = k_list[np.where(pf == np.max(pf))[0][0]]
        opt_criteria = list(self.criteria_dict.keys())[np.where(pf == np.max(pf))[1][0]]

        skb = SelectKBest(self.criteria_dict[opt_criteria], k=opt_k)
        skb.fit(self.X_train, self.y_train.astype(np.float).ravel())
        self.skb = skb

        return opt_k, opt_criteria, opt_C, roc_auc
    
    def transform(self, X):
        select_X = self.skb.transform(X)
        return select_X
    
    def extract_genes(self, colnames, n_genes=500): # takes list of column names corresponding to the original input matrix
        gene_set = colnames[self.skb.scores_.argsort()[-500:]]

        # separate genes which are connected together with a ";"
        sep_geneset = []
        for i in range(len(gene_set)):
            split_list = gene_set[i].split(";") # element in list can contain multiple gene names (overlapping)
            for gene in split_list:
                if gene in sep_geneset:
                    continue
                else:
                    sep_geneset.append(gene)

        return sep_geneset
