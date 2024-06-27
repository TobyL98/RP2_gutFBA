####################
# species_choice.py
###################

import pandas as pd
from pathlib import Path


def top_species(df):

    '''Picks out the species in the genus that is on average most abundant'''
    df['average_abundance'] = df.iloc[:, 1:].mean(axis = 1)
    average_df = df.loc[:, ['species', 'average_abundance']]

    top_df = average_df.groupby('Genus').max()
    top_df = top_df.sort_values(by = ['average_abundance'], ascending = False)

    return average_df, top_df


# read in the healthy and stage_I_II stages
healthy_path = Path('Data/top_genushealthy.csv')
healthy_df = pd.read_csv(healthy_path, sep = ',', index_col= 'Genus')
healthy_df = healthy_df.drop(columns = ['Unnamed: 0'])

stage_I_II_path = Path('Data/top_genusCRC.csv')
stage_I_II_df = pd.read_csv(stage_I_II_path, sep = ",", index_col = 'Genus')
stage_I_II_df = stage_I_II_df.drop(columns = ['Unnamed: 0'])


averagegenus_healthy_df, topspecies_healthy_df = top_species(healthy_df)
topspecies_healthy_df.to_csv("model_picker/topspecies_healthy_df.csv")
averagegenus_healthy_df.to_csv("model_picker/averagegenus_healthy_df.csv")

averagegenus_stageI_II_df, topspecies_stageI_II_df = top_species(stage_I_II_df)
topspecies_stageI_II_df.to_csv("model_picker/topspecies_stageI_II_df.csv")
averagegenus_stageI_II_df.to_csv("model_picker/averagegenus_stageI_II_df.csv")
