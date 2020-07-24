import pandas as pd

def DemoData(path_to_demographic, patid_list):
    demo_sourcedf = pd.read_csv(path_to_demographic)
    demo_df = pd.DataFrame(index=patid_list, columns=["Age", "Gender 1", "Gender 2", "Asian", "Black", "Islander", "White"])
    
    # for i in range(len(patid_list)):
    #     if type(patid_list[i]) == str:
    #         patid_list[i] = int(patid_list[i][-4:])
    #     elif type(patid_list[i]) == int:
    #         continue
    #     else:
    #         print("Patient ID is unexpected: neither integer nor string.")
    
    for patid in patid_list:
        int_patid = int(patid[-4:])
        record = demo_sourcedf[demo_sourcedf["PATNO"] == int_patid]
        
        # Get age
        last_update = int(record["LAST_UPDATE"].values[0][:4])
        birth_year = int(record["BIRTHDT"].values[0]) 
        age = last_update - birth_year
        
        # Get gender
        gender_1 = (record["GENDER"] == 1).sum()
        gender_2 = (record["GENDER"] == 2).sum()

        # Get race
        asian = record["RAASIAN"].values[0]
        black = record["RABLACK"].values[0]
        islander = record["RAHAWOPI"].values[0]
        white = record["RAWHITE"].values[0]

        # Put it together
        row_dict = {"Age" : age, "Gender 1" : gender_1, "Gender 2" : gender_2, "Asian" : asian, "Black" : black, "Islander" : islander, "White" : white}
        for key in row_dict:
            demo_df.loc[patid][key] = row_dict[key]
        
    return demo_df
