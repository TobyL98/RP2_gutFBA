##################
# flux_getter.py
#################

'''Gets a specific secretion flux from each individual sample'''

import pandas as pd
from pathlib import Path

def get_flux(sec_flux, samp_dir):

    '''Gets a specific secretion flux
    from each sample'''
    reaction_sec_dict = {}
    for sample in samp_dir.rglob('*'):
        sample_stem = sample.stem
        if sample_stem.startswith("sec"):
            sample_number = sample_stem.split('_')[3]

            sample_sec_df = pd.read_csv(sample)
            sample_reaction_sec_df = sample_sec_df.loc[sample_sec_df['reaction'] == sec_flux, 'flux']
            reaction_sec_num = sample_reaction_sec_df.values[0]
            
            reaction_sec_dict[sample_number] = reaction_sec_num
    all_reaction_sec_df = pd.DataFrame.from_dict(reaction_sec_dict,
                                               orient= 'index',
                                               columns= ['Secretion Flux'])
    print(all_reaction_sec_df)
        


secretion_flux = "EX_ac_medium"
samples_dir = Path('variation_healthy')

get_flux(secretion_flux, samples_dir)


        