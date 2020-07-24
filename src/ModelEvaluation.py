import numpy as np
from sklearn import metrics
from sklearn.metrics import auc
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
from numpy import interp
from itertools import cycle
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import StratifiedKFold


class ROCMulti:
    """
    This class was built to handle processes around model evaluation for PD progression subtype classification.
    Arguments: 
        model, X, y, n_classes (default 3)
    
    Methods: 
        calculate macro-averaged tpr/fpr and ROC-AUC, plot 
        calculate micro-averaged tpr/fpr and ROC-AUC, plot 
        calculate OVR tpr/fpr and ROC-AUC for class X, plot
        plot all ROCs mentioned above
    
    Usage:
        rocmulti = ROCMulti(model, X, y)
        macro_auc = rocmulti.macro_auc()
        micro_auc = rocmulti.micro_auc()
        ovr_auc = rocmulti.ovr_auc(class_label=1)
        rocmulti.plot_all()
    """

    def __init__(self, model, X, y, n_classes=3):
        self.input_data = X
        self.target_data = y
        self.n_classes = n_classes

        # get continuous value for prediction based on model –– either decision function or probability
        try:
            self.y_score = model.decision_function(X)
        except:
            self.y_score = model.predict_proba(X)
        
        # next create a one-hot coded vector of the target array
        self.y_onehot = OneHotEncoder().fit_transform(y.reshape(-1,1)).toarray()


    def macro_auc(self):
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(self.n_classes):
            fpr[i], tpr[i], _ = roc_curve(self.y_onehot[:, i], self.y_score[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])
        all_fpr = np.unique(np.concatenate([fpr[i] for i in range(self.n_classes)]))
        # interpolate all ROC curves at this points
        mean_tpr = np.zeros_like(all_fpr)
        for i in range(self.n_classes):
            mean_tpr += interp(all_fpr, fpr[i], tpr[i])
        # average it and compute AUC
        mean_tpr /= self.n_classes
        roc_auc = auc(all_fpr, mean_tpr)
        return roc_auc
    
    def micro_auc(self):
        fpr, tpr, _ = roc_curve(self.y_onehot.ravel(), self.y_score.ravel())
        roc_auc = auc(fpr, tpr)
        return roc_auc
    
    def ovr_auc(self, class_label=None): # Either returns dict of all classes or single value of specified class
        if class_label is None:
            fpr = dict()
            tpr = dict()
            roc_auc = dict()
            for i in range(self.n_classes):
                fpr[i], tpr[i], _ = roc_curve(self.y_onehot[:, i], self.y_score[:, i])
                roc_auc[i] = auc(fpr[i], tpr[i])
        else:
            fpr, tpr, _ = roc_curve(self.y_onehot[:, class_label], self.y_score[:, class_label])
            roc_auc = auc(fpr, tpr)
        return roc_auc
    
    def plot_all(self):
        # Compute ROC curve and ROC area for each class
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(self.n_classes):
            fpr[i], tpr[i], _ = roc_curve(self.y_onehot[:, i], self.y_score[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])
        
        # compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(self.y_onehot.ravel(), self.y_score.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

        all_fpr = np.unique(np.concatenate([fpr[i] for i in range(self.n_classes)]))
        
        # interpolate all ROC curves at this points
        mean_tpr = np.zeros_like(all_fpr)
        for i in range(self.n_classes):
            mean_tpr += interp(all_fpr, fpr[i], tpr[i])
        
        # average it and compute AUC
        mean_tpr /= self.n_classes

        fpr["macro"] = all_fpr
        tpr["macro"] = mean_tpr
        roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])
        
        # plot all ROC curves
        plt.figure(figsize=(8,8))
        plt.plot(fpr["micro"], tpr["micro"],
                label="micro-average ROC curve (area = {0:0.2f})"
                    "".format(roc_auc["micro"]),
                color="deeppink", linestyle=":", linewidth=4)

        plt.plot(fpr["macro"], tpr["macro"],
                label="macro-average ROC curve (area = {0:0.2f})"
                    "".format(roc_auc["macro"]),
                color="navy", linestyle=":", linewidth=4)

        colors = cycle(["aqua", "darkorange", "cornflowerblue"])
        for i, color in zip(range(self.n_classes), colors):
            plt.plot(fpr[i], tpr[i], color=color, lw=2,
                    label="ROC curve of class {0} (area = {1:0.2f})"
                    "".format(i+1, roc_auc[i]))
        
        plt.plot([0, 1], [0, 1], "k--", lw=2)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("ROC of PD subtype prediction")
        plt.legend(loc="lower right")



#####

def roc_binary(model, X, y): # , nr=5
    plt.figure(0, figsize=(8,8))
    try:
        y_decf = model.decision_function(X)
    except:
        y_decf = model.predict_proba(X)
    
    fpr, tpr, _ = metrics.roc_curve(y, y_decf)
    auc = metrics.auc(fpr, tpr)
    plt.plot(
        fpr, tpr, 
        label="ROC curve (area = {0:0.2f})".format(auc),
        color="deeppink", linestyle=":", linewidth=4
        )

    plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--",)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC of PD binary prediction")
    plt.legend(loc="lower right")

    # aucs = np.zeros(nr)

    # for r in range(nr):
    #     skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=r*8800+1)
    #     clf = model(params, class_weight="balanced", max_iter=3000, n_jobs=6)
    #     y_true = np.array([])
    #     y_decf = np.array([])
    #     for train_index, test_index in skf.split(X, y):
    #         X_train, X_test = X[train_index,:], X[test_index,:]
    #         y_train, y_test = y[train_index], y[test_index]
    #         clf.fit(X_train, y_train)
    #         decf = clf.decision_function(X_test)
    #         y_true = np.concatenate((y_true, y_test))
    #         y_decf = np.concatenate((y_decf, decf))
    #     fpr, tpr, _ = metrics.roc_curve(y_true, y_decf)
    #     aucs[r] = metrics.auc(fpr, tpr)
    #     plt.plot(fpr, tpr, lw=lw)

    # plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    # plt.xlabel("False Positive Rate")
    # plt.ylabel("True Positive Rate")
    # plt.title("Logistic Regression with C=1e+5 (mean AUC = {:.3f})".format(np.mean(aucs)), size=16)


#####
