#####################
# shannon_diversity
#####################

'''For the CRC and healthy genus and species data.
Calculates the Shannon diversity and counts the number of species'''

import pandas as pd
import numpy as np
from pathlib import Path
import math

### FUNCTIONs
def shannon(col):

    '''Takes each column with bacterial abundnace as an input
      and creates a new column calculating Shannon index.
      Will also count the number of species >0 to calculate richness'''
    col = col.loc[col != 0.0]
    proportion_col = col / sum(col)
    shannon_col = proportion_col * np.log(proportion_col)
    shannon_index = -sum(shannon_col)
    richness = len(col)

    
    return(shannon_index, richness)

def apply_shannon(fp):

    '''Takes in the filepath and calculates
    the shannon index from the dataframe given'''
    df = pd.read_csv(fp, sep= ',')
    shannon_df = (df
                  .iloc[:, 2:]
                  .apply(shannon, axis= 'index')
                  .transpose()
                  .rename(columns= {0: 'shannon_index',
                                    1: 'richness'}))
    return(shannon_df)

# run for genus healthy
genus_healthy_path = Path('Outputs/healthy_df_out.csv')
genus_healthy_shannon_df = apply_shannon(genus_healthy_path)
genus_healthy_output_path = Path('graphs/diversity/genus_healthy.csv')
genus_healthy_shannon_df.to_csv(genus_healthy_output_path)

# run for genus CRC
genus_CRC_path = Path('Outputs/Stage_I_II_df_out.csv')
genus_CRC_shannon_df = apply_shannon(genus_CRC_path)
genus_CRC_output_path = Path('graphs/diversity/genus_CRC.csv')
genus_CRC_shannon_df.to_csv(genus_CRC_output_path)

# run for species healthy
species_healthy_path = Path('Data/species_healthy.csv')
species_healthy_shannon_df = apply_shannon(species_healthy_path)
species_healthy_output_path = Path('graphs/diversity/species_healthy.csv')
species_healthy_shannon_df.to_csv(species_healthy_output_path)

# run for species CRC
species_CRC_path = Path('Data/species_CRC.csv')
species_CRC_shannon_df = apply_shannon(species_CRC_path)
species_CRC_output_path = Path('graphs/diversity/species_CRC.csv')
species_CRC_shannon_df.to_csv(species_CRC_output_path)

# calculating shannon and richness for healthy at 95%
healthy95_path = Path('Outputs/average/healthy_df_out_ave95.csv')
healthy95_df = pd.read_csv(healthy95_path, sep= ',')
healthy95_df = pd.DataFrame(healthy95_df['average_abundance'])
healthy95_shannon_df = (healthy95_df
                      .apply(shannon)
                      .transpose()
                      .rename(columns= {0: 'shannon_index',
                                    1: 'richness'}))
healthy95_output_path = Path('graphs/diversity/healthy95.csv')
healthy95_shannon_df.to_csv(healthy95_output_path)

# calculating shannon and richness for CRC at 95%
crc95_path = Path('Outputs/average/StageI_II_df_out_ave95.csv')
crc95_df = pd.read_csv(crc95_path, sep= ',')
crc95_df = pd.DataFrame(crc95_df['average_abundance'])
crc95_shannon_df = (crc95_df
                      .apply(shannon)
                      .transpose()
                      .rename(columns= {0: 'shannon_index',
                                    1: 'richness'}))
crc95_output_path = Path('graphs/diversity/crc95.csv')
crc95_shannon_df.to_csv(crc95_output_path)






