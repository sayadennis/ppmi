import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn import metrics


# Logistic Regression: Learn hyperparameters by cross-validation
def opt_hyp_lrm(X, y, params_list, multi=True):
    if multi:
        multi_class = "ovr"
        scoring = "roc_auc_ovr_weighted"
        refit = "roc_auc_ovr_weighted"
    else:
        multi_class = "auto"
        scoring = "roc_auc"
        refit = "roc_auc"
    
    lrCV = LogisticRegressionCV(class_weight="balanced", n_jobs=8, max_iter=3000,
        Cs=params_list, multi_class=multi_class, scoring=scoring, cv=StratifiedKFold(n_splits=5, random_state=1, shuffle=True), refit=refit
    )
    lrCV.fit(X,y)
    opt_C = lrCV.C_[0]
    opt_mean_score = np.mean(lrCV.scores_[1][:,np.where(params_list == opt_C)])

    return opt_C, opt_mean_score


# SVM: Learn hyperparameters by cross-validation
def opt_hyp_svm(X, y, params_list, multi=True):
    if multi:
        scoring = "roc_auc_ovr_weighted"
        refit = "roc_auc_ovr_weighted"
    else:
        scoring = "roc_auc"
        refit = "roc_auc"

    gsCV = GridSearchCV(
        SVC(class_weight="balanced", max_iter=3000, random_state=0, probability=True, decision_function_shape="ovr"), n_jobs=8, 
        param_grid=params_list, scoring=scoring, cv=StratifiedKFold(n_splits=5, random_state=1, shuffle=True), refit=refit
    )
    gsCV.fit(X,y)
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

