###############
# flux_all_picker.py
###############

'''Picks all fluxes with an average greater 1e-5
and calculates the median and IQR for the samples of 10'''

import pandas as pd
import numpy as np
from pathlib import Path
from flux_picker import get_flux, flux_finder

def flux_all_picker(sec_flux_df):

    '''Picks out all fluxes < -1e-5 that can be 
    compared against the individual analysis.'''

    sec_flux_df = sec_flux_df.loc[sec_flux_df['flux'] < -1e-5]
    print(sec_flux_df.head())
    reaction_avoid = ["community_biomass", "EX_h2_medium", "EX_h_medium", "EX_co2_medium"]
    sec_flux_df = sec_flux_df.loc[~sec_flux_df['reaction'].isin(reaction_avoid)]
    flux_reactions = sec_flux_df['reaction'].to_list()

    return flux_reactions, sec_flux_df

def main():
    # run with healthy data
    healthy_average_flux_path = Path('Western_healthy/sec_flux_Western_healthy.csv')
    healthy_average_flux_df = pd.read_csv(healthy_average_flux_path, sep= ',')
    flux_reactions, flux_df = flux_all_picker(healthy_average_flux_df)

    healthy_samples_dir = Path("variation_healthy")
    print('Healthy')
    flux_finder(flux_reactions, flux_df, healthy_samples_dir, data= "healthy")

    # run with CRC data
    CRC_average_flux_path = Path('Western_CRC/sec_flux_Western_CRC.csv')
    CRC_average_flux_df = pd.read_csv(CRC_average_flux_path, sep = ',')
    # use random flux reactions to get CRC random df
    CRC_flux_df = CRC_average_flux_df.loc[CRC_average_flux_df['reaction'].isin(flux_reactions)
                                                 , :]
    CRC_samples_dir = Path('variation_CRC')
    print('CRC')
    flux_finder(flux_reactions, CRC_flux_df, CRC_samples_dir, data= 'crc')


if __name__ == "__main__":
    main()
