
def applyThres(table, thres=0.05):
    # table: Pandas DataFrame
    # thres: float
    bool_col = [] # list of boolean indicating the columns to keep at True

    for colname in table.columns:
        ct = 0
        for i in range(len(table[colname])):
            if table[colname][i] > 0:
                ct += 1
            else:
                continue
        rate = ct/len(table[colname])
        if ((rate < thres) | (rate > (1-thres))):
            bool_col.append(False)
        else:
            bool_col.append(True)
    
    return table.iloc[:,bool_col]
