###########################
# MES.py
###########################

'''Will calculate the metabolic exchange
score. This is based on the results 
from the potential_metabolite_exchanges
in python.'''

import pandas as pd
from pathlib import Path
import re

def counter(species):
    '''Counts the number of bacterial
    species in produced by or consumed by'''
    
    species_count = len(re.split(r',', species))
    return species_count


def MES(path):

    '''Function reads in the dataframe, calculates the metabolic
    exchange score and adds it as a column to the dataframe'''

    d_type = {'metabolite_id': str,
          'metabolite_name': str,
          'cross_feeding': bool,
          'produced_by': str,
          'consumed_by': str}
    use_cols = ['metabolite_id', 'metabolite_name', 'cross_feeding', 'produced_by', 'consumed_by']
    df = pd.read_csv(path, usecols= use_cols, dtype = d_type)


    df = df.loc[df['cross_feeding'] == True]
    df.loc[:, ['produced_count']] = df['produced_by'].apply(counter)
    df.loc[:, ['consumed_count']] = df['consumed_by'].apply(counter)

    # calculating MES
    df.loc[:, ['MES']] = 2 * ((df['produced_count'] * df['consumed_count']) / (df['produced_count'] + df['consumed_count']))
    df = df.sort_values(by = ['MES'], ascending = False)
    
    return(df)

def MES_difference(df1, df2):

    '''Calculates the difference in
    MES between two groups'''

    # merge the two groups
    MES_merged_df = df1.merge(df2, on = 'metabolite_id', how = 'outer')
    MES_merged_df['MES_diff'] = MES_merged_df['MES_health'] - MES_merged_df['MES_CRC']
    MES_merged_df = MES_merged_df.sort_values(by =['MES_diff'], ascending = True)
    print(MES_merged_df.head())

    MES_merged_df.to_csv("MES_final_results.csv")



def main():

    '''Runs the full function. Inputs are
    the metabolic exchnage results from PyCoMo
    for two species.'''
    healthy_met_exchange_path = Path("healthy/mes_results.csv")
    healthy_MES_df = MES(healthy_met_exchange_path)

    CRC_met_exchnage_path = Path("StageI_II/mes_results.csv")
    CRC_MES_df = MES(CRC_met_exchnage_path)

    # rename df columns prior to the merge
    healthy_rename_dict = {"cross_feeding": "cross_feeeding_health",
                           "produced_by": "produced_by_health",
                           "consumed_by": "consumed_by_health",
                           "produced_count": "produced_count_health",
                           "consumed_count": "consumed_count_health",
                           "MES": "MES_health"}
    healthy_MES_df = healthy_MES_df.rename(columns= healthy_rename_dict)

    CRC_rename_dict = {"cross_feeding": "cross_feeeding_CRC",
                           "produced_by": "produced_by_CRC",
                           "consumed_by": "consumed_by_CRC",
                           "produced_count": "produced_count_CRC",
                           "consumed_count": "consumed_count_CRC",
                           "MES": "MES_CRC"}
    CRC_MES_df = CRC_MES_df.rename(columns = CRC_rename_dict)
    MES_difference(healthy_MES_df, CRC_MES_df)

main()