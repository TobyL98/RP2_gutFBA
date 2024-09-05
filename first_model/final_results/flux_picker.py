##################
# flux_picker
##################

'''Picks 3/4 fluxes for further analysis'''

import pandas as pd
import numpy as np
from pathlib import Path

def get_flux(sec_reaction, samp_dir, reaction_secretion_dict):

    '''Gets a specific secretion flux
    from each sample'''

    for sample in samp_dir.rglob('*'):
        sample_stem = sample.stem
        # find secretion flux files
        if sample_stem.startswith("sec"):
            sample_number = sample_stem.split('_')[3]

            # gets the flux value and adds to dictionary
            sample_sec_df = pd.read_csv(sample)
            sample_reaction_sec_df = sample_sec_df.loc[sample_sec_df['reaction'] == sec_reaction, 'flux']
            if sample_reaction_sec_df.shape[0] != 0:
                reaction_sec_flux = sample_reaction_sec_df.values[0]
                reaction_secretion_dict[sample_number] = reaction_sec_flux
    

    return reaction_secretion_dict

def flux_picker(sec_flux_df):

    '''Picks out four random secretion fluxes that can be 
    compared against the individual analysis.
    Then uses function getr_flux to create dataframes with
    the secretion flux from each individual analysis'''

    sec_flux_df = sec_flux_df.loc[sec_flux_df['flux'] < -1e-5]


    # pick 4 random fluxes to sample
    np.random.seed(42)
    overall_fluxes_num = sec_flux_df.shape[0]
    fluxes_range = range(1, overall_fluxes_num)
    random_fluxes = np.random.choice(fluxes_range, size= 4, replace= False)
    rand_flux_df = sec_flux_df.iloc[random_fluxes, :]
    random_flux_reactions = rand_flux_df['reaction'].to_list()

    return random_flux_reactions, rand_flux_df


def flux_finder(random_flux_reactions, rand_flux_df, samples_dir, data):

    '''For each reaction in the random flux dataframe.
    Will get the flux for every sample from healthy and CRC
    variation and add it to a healthy or CRC dataframe
    '''
    reaction_sec_dict = {}

    # goes through four random reactions and access
    # the flux values from the individual samples for each
    for reaction in random_flux_reactions:
        reaction_sec_dict = {}
        reac_flux = (rand_flux_df
                     .loc[rand_flux_df['reaction'] == reaction, 'flux']
                     .values[0])
        reaction_sec_dict['average'] = reac_flux
        reaction_sec_dict = get_flux(reaction, samples_dir, reaction_sec_dict)
        all_reaction_sec_df = (pd.DataFrame
                               .from_dict(reaction_sec_dict,
                                               orient= 'index',
                                               columns= ['flux']))
        all_reaction_sec_df.index.name = 'sample'
        # save reaction_df for all samples run
        if data == 'healthy':
            csv_path = samples_dir / 'variation_{0}_healthy.csv'.format(reaction)
        else:
            csv_path = samples_dir / 'variation_{0}_crc.csv'.format(reaction)
        print(all_reaction_sec_df)
        all_reaction_sec_df.to_csv(csv_path)


    # run get_flux with CRC data
def main():
    healthy_average_flux_path = Path('Western_healthy/sec_flux_Western_healthy.csv')
    healthy_average_flux_df = pd.read_csv(healthy_average_flux_path, sep= ',')
    random_flux_reactions, healthy_random_flux_df = flux_picker(healthy_average_flux_df)

    # run with healthy data
    health_samples_dir = Path('variation_healthy')
    print('Healthy')
    flux_finder(random_flux_reactions, 
                healthy_random_flux_df, 
                health_samples_dir,
                data = "healthy")

    # run with CRC data
    CRC_average_flux_path = Path('Western_CRC/sec_flux_Western_CRC.csv')
    CRC_average_flux_df = pd.read_csv(CRC_average_flux_path, sep = ',')
    # use random flux reactions to get CRC random df
    CRC_random_flux_df = CRC_average_flux_df.loc[CRC_average_flux_df['reaction'].isin(random_flux_reactions)
                                                 , :]
    CRC_samples_dir = Path('variation_CRC')
    print('CRC')
    flux_finder(random_flux_reactions, CRC_random_flux_df, CRC_samples_dir, data= "CRC")

if __name__ == "__main__":
    main()



