#################
# Abundance calculator
#################

# Code will calculate the abundances of
# data from Yachida et al.,
# See the data in data folder
# can then be used in FBA analysis

import pandas as pd
from pathlib import Path
import sys

def read_sort():

    '''Function reads in the two initial data files.
    Ones with the groups the patients are in (e.g., helathy, Stage I and II cancer)
    The other with the levels of bacteria in the gut from metagenomics.
    Finally, it sorts the two data files'''

    # read in the data
    samples_group_path = Path('Data/Yachida_data_groups.csv')
    use_cols = ["Subject_ID", "Group"]
    d_type = {"Subject_ID": 'category',
              "Group": 'category'}
    sample_groups_df = pd.read_csv(samples_group_path, sep =',', usecols= use_cols, dtype= d_type)

    sample_levels_path = Path('Data/Yachida_data_levels.csv')
    sample_levels_df = pd.read_csv(sample_levels_path, sep = ',')

    #organising sample_levels_df
    sample_levels_df = sample_levels_df.rename(columns= {"Unnamed: 0": "Species"})

    # remove any viruses or phage as not wanted
    # also removes some known protozoa and other things not needed in dataset
    sample_levels_df["virus"] = sample_levels_df["Species"].str.contains(
        r"([Vv]irus)|([Pp]hage)|([Gg]iardia)|([Ss]accharomyces)|([Cc]andidate)", 
        regex= True
        )
    sample_levels_df = sample_levels_df.loc[sample_levels_df['virus'] == False]
    sample_levels_df = sample_levels_df.drop(columns= 'virus')



    return sample_groups_df, sample_levels_df


def merge_filter(group_df, level_df):

    '''Function merges the group and level data on Subject_ID.
    Then filters by 'Healthy' and 'Stage_I_II'''''

    # pivot so headers are species
    # so it can join with sample_groups_df
    level_df = level_df.set_index('Species')
    level_df = level_df.transpose()

    # rename as now Subject_ID rather than species
    level_df = level_df.rename(columns= {"Species": "Subject_ID"})

    # merge the two tables
    merged_df = pd.merge(group_df, level_df, left_on = "Subject_ID", right_index = True, how = 'inner')

    # Filter to get two tables
    # One for healthy and one for Stage_I_II
    healthy_df = merged_df.loc[merged_df['Group'] == "Healthy"]
    Stage_I_II_df = merged_df.loc[merged_df['Group'] == "Stage_I_II"]

    # drop Group column as no longer needed
    healthy_df = healthy_df.drop(columns = ['Group'])
    Stage_I_II_df = Stage_I_II_df.drop(columns = ['Group'])

    return healthy_df, Stage_I_II_df

def species_to_genus(species):
    genus = species.split(" ")[0]
    return(genus)

def group_genus(df, ID):
    '''Function creates a genus column from
    the species column. Then groups by the genus column.'''

    df = df.set_index('Subject_ID')
    df = df.transpose()

    # reset index so species names can be edited
    df = df.reset_index()
    df = df.rename(columns = {'index': 'species'})

    df['Genus'] = df['species'].apply(species_to_genus)
    df.to_csv("top_genus{0}.csv".format(ID))

    # group by new genus column and sum each species abundance
    df = df.groupby(['Genus']).sum()

    return(df)

def main():

    '''Runs the main code'''

    # obtains initial data and sorts it
    # one df for groups the patient is in (e.g., helathy, stage 1 cancer)
    # one df for levels of bacteria
    groups_df, levels_df = read_sort()

    # merges the two tables to combine group data and lecvels of bacteria
    # filters to create and healthy and Stage I and II CRC group
    healthy_df, Stage_I_II_df = merge_filter(groups_df, levels_df)

    # create the genus column in healthy df
    healthy_df = group_genus(healthy_df, "healthy")
    print(healthy_df.head())
    healthy_df.to_csv("Outputs/healthy_df_out.csv")

    # create the genus column in Stage_I_II_df
    Stage_I_II_df = group_genus(Stage_I_II_df, "CRC")
    print(Stage_I_II_df.head())
    Stage_I_II_df.to_csv("Outputs/Stage_I_II_df_out.csv")

if __name__ == "__main__":
    sys.exit(main())
