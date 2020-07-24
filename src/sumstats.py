import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_hist(X):
    # Histogram of data (both for columns and rows) –– THIS IS USES PRE-NORMALIZED DATA TO GET ACTUAL COUNTS
    fig, axs = plt.subplots(1, 2, figsize=(14,4))
    axs[0].hist(np.sum(X, axis=0), bins=20)
    axs[0].set_title("Number of patients for each variant", size=14)
    axs[1].hist(np.sum(X, axis=1), bins=20)
    axs[1].set_title("Number of variants in a patient", size=14)


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

