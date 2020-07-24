import pandas as pd

class AlignInputTarget:
    """
    This class aligns the input and target dataframes of the PD datasets based on patient ID.
    The function is written to align variant count data as input, the indices of which are PPMI_SI_<num> (string),
    and PD progression subtype as target, the indices of which are simply <num> (integer).
    Each method returns the appropriately reduced input dataframe and target dataframe.

    Example of how to use this:
    aligner = AlignInputTarget(X, y)
    aligned_input = aligner.align_input()
    aligned_target = aligner.align_target()
    """
    
    def __init__(self, input_df, target_df):
        self.inputindex = input_df.index
        self.targetindex = target_df.index
        self.input_df = input_df
        self.target_df = target_df
        self.ol_patid = []
        for i in range(len(self.targetindex)):
            for j in range(len(self.inputindex)):
                if str(self.targetindex[i]) in self.inputindex[j]:
                    self.ol_patid.append(self.targetindex[i])
                else:
                    continue
    
    def align_input(self): # keeps the input's index order 
        bool_input = []
        for patno in self.inputindex:
            if int(patno[-4:]) in self.ol_patid:
                bool_input.append(True)
            else:
                bool_input.append(False)
        new_input = self.input_df[bool_input]
        return new_input
    
    def align_target(self): # changes the target's index order according to the input's index order
        new_target = pd.DataFrame(columns=["CLUSTER_IDX"])
        for patno in self.inputindex:
            if int(patno[-4:]) in self.ol_patid:
                new_target = new_target.append(
                    pd.DataFrame(self.target_df.loc[int(patno[-4:])].values, index=[int(patno[-4:])], columns=["CLUSTER_IDX"])
                    )
            else:
                continue
            
        # for patno in self.targetindex:
        #     if patno in self.ol_patid:
        #         self.bool_target.append(True)
        #     else:
        #         self.bool_target.append(False)
        # new_target = self.target_df[self.bool_target]
        return new_target
    