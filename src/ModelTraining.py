import numpy as np
from sklearn import metrics

from sklearn.linear_model import LogisticRegressionCV
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV

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
        opt_mean_score = np.mean(lrCV.scores_[1][:,np.where(param_list == opt_C)])

        return opt_C, opt_mean_score
    
    def tune_svm(self, X, y, param_list=None):
        # set default parameter grid
        if param_list == None:
            param_list = self.svm_params

        # specify scoring and refit criteria respectively for multiclass and binary
        if self.multi:
            scoring = "roc_auc_ovr_weighted"
            refit = "roc_auc_ovr_weighted"
        else:
            scoring = "roc_auc"
            refit = "roc_auc"

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



class TuneSelectKBest:
    """
    This class handles the tuning of the parameter K as well as selection criteria for SelectKBest function.
    """

    def __init__(self):
        self.criteria_dict = {
            "f_classif" : f_classif, 
            "mutual_info_classif" : mutual_info_classif, 
            "chi2" : chi2, 
        }
        self.y = y
    
    def tune_selectkbest(self, k_list, hypertune_func, multi=True, eval="micro"): # eval can be "micro" or "macro"
        best_k = 0
        best_criteria = 0
        return best_k, best_criteria
    
    def transform(self, X_test):
        transformed_X_test = 0
        return transformed_X_test
    
    def extract_genes(self):
        gene_set = []
        return gene_set

