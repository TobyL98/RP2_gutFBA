#####################
# flux_change.py
#####################

'''Calculates whether the flux has changed from
an uptake to a secretion flux or vice versa
between the CRC and healthy dataset'''

import pandas as pd
from pathlib import Path

def sec_flux_diff():
    flux_healthy_path = Path('StageI_II/sec_flux.csv')
    d_type = {'flux': float,
              'reaction': 'category',
              'metabolite': 'category'}
    use_cols = ["flux", "reaction", "metabolite"]
    flux_healthy_df = pd.read_csv(flux_healthy_path, sep = ",", dtype = d_type, usecols = use_cols)


    flux_CRC_path = Path('healthy/sec_flux.csv')
    flux_CRC_df = pd.read_csv(flux_CRC_path, sep = ",", dtype = d_type, usecols = use_cols)

    merged_df = pd.merge(flux_healthy_df, flux_CRC_df, on = ["reaction", "metabolite"], how = 'outer')
    merged_df = merged_df.rename(columns= {"flux_x": "flux_healthy",
                                       "flux_y": "flux_CRC"})
    merged_df = merged_df.loc[:,['metabolite', 'reaction', 'flux_healthy', "flux_CRC"]]

    # find where the NaN values are in any row, should show where changes are
    is_na_df = merged_df.isna().any(axis = 1)
    secCRC_to_uptakeHealthy_df = merged_df[is_na_df]
    print(secCRC_to_uptakeHealthy_df)




sec_flux_diff()