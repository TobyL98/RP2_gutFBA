##################
# flux_picker
##################

'''Picks 3/4 fluxes for further analysis'''

import pandas as pd
import numpy as np

sec_flux_df = pd.read_csv('sec_flux_Western_healthy.csv', sep= ',')

def flux_picker():

    '''Picks out four random secretion fluxes that can be 
    compared against the individual analysis'''

    sec_flux_df = sec_flux_df.loc[sec_flux_df['flux'] < -1e-5]
    # pick 4 random fluxes to sample
    np.random.seed(42)
    overall_fluxes_num = sec_flux_df.shape[0]
    fluxes_range = range(1, overall_fluxes_num)
    random_fluxes = np.random.choice(fluxes_range, size= 4, replace= False)
    rand_flux_df = sec_flux_df.iloc[random_fluxes, :]
    print(rand_flux_df)

