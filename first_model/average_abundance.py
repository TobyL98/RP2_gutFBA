##################
# average_abundance.py
##################

import pandas as pd
from pathlib import Path


'''Calculates the average abundance of bacterial species
for each participant. Then calculates the cumalative abundance
and picks out the number of species at a certain cumalative abundance
e.g., 0.95 (95%) or 0.99 (99%)'''

# read in files
healthy_path = Path("Outputs/healthy_df_out.csv")
healthy_df = pd.read_csv(healthy_path, sep = ',')

Stage_I_II_path = Path("Outputs/Stage_I_II_df_out.csv")
Stage_I_II_df = pd.read_csv(Stage_I_II_path, sep = ',')

def average_abundance(df, cut_off):

    '''Function will calculate the average
    abudance of the bacterial species.
    Will also create a column for cumalative abundance'''

    df['average_abundance'] = df.iloc[:, 2:].mean(axis = 1)
    df['average_abundance'] = df['average_abundance'] / sum(df['average_abundance'])
    
    # remove individual columns
    df = df.loc[:, ['Genus', 'species', 'average_abundance']]
    
    # cumalative abundance
    df = df.sort_values(by = ['average_abundance'], ascending = False)
    df['cumalative_abundance'] = df['average_abundance'].cumsum()
    
    df.to_csv('checks.csv')
    df = df.loc[df['cumalative_abundance'] <= cut_off]

    return df


# main
cut_off = 0.35

healthy_average_df = average_abundance(healthy_df, cut_off)
StageI_II_average_df = average_abundance(Stage_I_II_df, cut_off)

healthy_average_df.to_csv("Outputs/average/healthy_df_out_ave_35.csv")
StageI_II_average_df.to_csv("Outputs/average/StageI_II_df_out_ave_35.csv")




