####################
# species_choice.py
###################

import pandas as pd
from pathlib import Path


def top_species(df):

    '''Picks out the species in the genus that is on average most abundant'''
    df['average_abundance'] = df.iloc[:, 1:].mean(axis = 1)
    df = df.loc[:, ['species', 'average_abundance']]

    # removing those that are unclassified 
    # and there is more than one species in genus
    # as won't be able to find unclassified
    # but keep unclassified if only one
    df['unclassified'] = df["species"].str.contains(r"[Uu]nclassified", regex = True)
    df['not_unique'] = df.index.duplicated(keep = False)
    df = df[(df['unclassified'] == False) | (df['not_unique'] == False)]
    df = df.drop(columns = ['unclassified', 'not_unique'])
    
    
    # picking species that is most abundant in the genus
    df = df.reset_index()
    df = df.sort_values(by = ['average_abundance'], ascending = False)
    top_df = df.drop_duplicates(subset = ['Genus'], keep = 'first')

    return top_df


# read in the healthy and stage_I_II stages
healthy_path = Path('Data/top_genushealthy.csv')
healthy_df = pd.read_csv(healthy_path, sep = ',', index_col= 'Genus')
healthy_df = healthy_df.drop(columns = ['Unnamed: 0'])

stage_I_II_path = Path('Data/top_genusCRC.csv')
stage_I_II_df = pd.read_csv(stage_I_II_path, sep = ",", index_col = 'Genus')
stage_I_II_df = stage_I_II_df.drop(columns = ['Unnamed: 0'])


topspecies_healthy_df = top_species(healthy_df)
topspecies_healthy_df.to_csv("model_picker/topspecies_healthydf.csv")

topspecies_stageI_II_df = top_species(stage_I_II_df)
topspecies_stageI_II_df.to_csv("model_picker/topspecies_stageI_IIdf.csv")
