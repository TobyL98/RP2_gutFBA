###################
# spec_diff.py
###################

'''Calculates the difference in the biomass flux
between species in healthy and CRC samples'''

import pandas as pd
from pathlib import Path

def perc_calc(df):
    
    '''Calculates the percentage contribution
    of ecah species to biomass. Outputs as new column'''

    total = sum(df['flux'])
    df['percentage_flux'] = (df['flux'] / total) * 100
    
    return df

def spec_diff(normal_df, cancer_df):

    '''Merges the healthy and CRC dataframe and
    calculates the difference in the percentage
    contribution to biomass flux for each species'''

    merged_df = pd.merge(normal_df, cancer_df, on = "Species_biomass", how = 'inner')
    merged_df = merged_df.rename(columns = {"flux_x": "flux_healthy",
                                            "flux_y": "flux_CRC",
                                            "percentage_flux_x": "percentage_flux_healthy",
                                            "percentage_flux_y": "percentage_flux_CRC"})
    merged_df['percentage_flux_diff'] = merged_df['percentage_flux_healthy'] - merged_df['percentage_flux_CRC']

    return merged_df

def spec_change(normal_df, cancer_df):

    '''worked out which species have changed from healthy to CRC
    i.e., which are not in the healthy or CRC sample'''

    normal_species = list(normal_df["Species_biomass"])
    CRC_species = list(cancer_df['Species_biomass'])
    
    print("Species in healthy but not CRC:")
    for species in normal_species:
        if species not in CRC_species:
            print(species)

    print("\nSpecies in CRC but not healthy:")
    for species in CRC_species:
        if species not in normal_species:
            print(species)



# read in samples
names = ['Species_biomass', 'flux']
healthy_path = Path("Western_healthy/spec_biomass_Western_healthy.csv")
healthy_df = pd.read_csv(healthy_path, 
                         sep = ",", 
                         names = names)

CRC_path = Path("Western_CRC/spec_biomass_Western_CRC.csv")
CRC_df = pd.read_csv(CRC_path, 
                     sep = ",",
                     names = names)

healthy_df = perc_calc(healthy_df)
CRC_df = perc_calc(CRC_df)

compare_df = spec_diff(healthy_df, CRC_df)
compare_df.to_csv("spec_flux_diff.csv")

spec_change(healthy_df, CRC_df)


