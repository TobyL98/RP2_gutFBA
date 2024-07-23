#######################
# flux_diff.py
#######################

'''Code will take in two dataframes. A dataframe for
uptake/secretion fluxes for both healthy and disease.
It can then merge the two dataframes and calculate what
are the biggest flux differences between healthy 
and disease.'''

import pandas as pd
from pathlib import Path


def flux_diff(healthy_df, CRC_df):

    merged_df = pd.merge(healthy_df, CRC_df, on = ["reaction", "metabolite"], how = 'inner')

    merged_df = merged_df.rename(columns= {"flux_x": "flux_healthy",
                                           "flux_y": "flux_CRC"})
    merged_df = merged_df.loc[:,['metabolite', 'reaction', 'flux_healthy', "flux_CRC"]]

    merged_df['flux_diff'] = merged_df['flux_healthy'] - merged_df['flux_CRC']
    merged_df['perc_flux_diff'] = (merged_df['flux_diff'] / 
                                   ((merged_df['flux_healthy'] + merged_df['flux_CRC']) / 2)) * 100
    merged_df = merged_df.sort_values(by = ['flux_diff'])

    print("Top 5 differences shifted towards healthy")
    print(merged_df.head())
    print("\nTop 5 differences shifted towards CRC")
    print(merged_df.tail())
    merged_df.to_csv('secFlux_diff.csv')


# read in data
flux_healthy_path = Path('Western_healthy/sec_flux_Western_healthy.csv')
d_type = {'flux': float,
          'reaction': 'category',
          'metabolite': 'category'}
use_cols = ["flux", "reaction", "metabolite"]
flux_healthy_df = pd.read_csv(flux_healthy_path, sep = ",", dtype = d_type, usecols = use_cols)


flux_CRC_path = Path('Western_CRC/sec_flux_Western_CRC.csv')
flux_CRC_df = pd.read_csv(flux_CRC_path, sep = ",", dtype = d_type, usecols = use_cols)

# run function
flux_diff(flux_healthy_df, flux_CRC_df)