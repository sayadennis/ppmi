import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_hist(X):
    # Histogram of data (both for columns and rows) –– THIS IS USES PRE-NORMALIZED DATA TO GET ACTUAL COUNTS
    fig, axs = plt.subplots(1, 2, figsize=(14,4))
    # number of patients with a variant per gene
    vals = []
    for i in range(X.shape[1]):
        ct = 0
        for j in range(len(X[:,i])):
            if X[j,i] > 0:
                ct += 1
            else:
                continue
        vals.append(ct)
    axs[0].hist(vals, bins=20)
    # axs[0].set_title("Number of patients per variant", size=14)
    axs[0].set_xlabel("Number of patients")
    axs[0].set_ylabel("Number of genes with variants")
    axs[1].hist(np.sum(X, axis=1), bins=20)
    # axs[1].set_title("Number of variants per patient", size=14)
    axs[1].set_xlabel("Number of variants")
    axs[1].set_ylabel("Number of patients")

def class_count(y):
    class_count = pd.DataFrame(0, index=["Subtype 1", "Subtype 2", "Subtype 3"], columns=["Counts"])
    for i in range(len(y)):
        if y[i] == 1:
            class_count.loc["Subtype 1"]["Counts"] += 1
        elif y[i] == 2:
            class_count.loc["Subtype 2"]["Counts"] += 1
        elif y[i] == 3:
            class_count.loc["Subtype 3"]["Counts"] += 1
        else:
            print("something went wrong")
    return class_count

