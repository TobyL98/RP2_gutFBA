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


def MES(df):

    '''Function calculates the metabolic
    exchange score'''

    df = df.loc[df['cross_feeding'] == True]
    df.loc[:, ['produced_count']] = df['produced_by'].apply(counter)
    df.loc[:, ['consumed_count']] = df['consumed_by'].apply(counter)

    # calculating MES
    df.loc[:, ['MES']] = 2 * ((df['produced_count'] * df['consumed_count']) / (df['produced_count'] + df['consumed_count']))
    df = df.sort_values(by = ['MES'], ascending = False)
    print(df.head())


met_exchange_path = Path("Toy/mes_results.csv")
d_type = {'metabolite_id': str,
          'metabolite_name': str,
          'cross_feeding': bool,
          'produced_by': str,
          'consumed_by': str}
use_cols = ['metabolite_id', 'metabolite_name', 'cross_feeding', 'produced_by', 'consumed_by']
met_exchange_df = pd.read_csv(met_exchange_path, usecols= use_cols, dtype = d_type)
MES(met_exchange_df)