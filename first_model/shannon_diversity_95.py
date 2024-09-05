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
def shannon_95(col):

    '''Takes each column with bacterial abundnace as an input
      and creates a new column calculating Shannon index.
      Will also count the number of species >0 to calculate richness'''
    col = col.loc[col != 0.0]
    col_df = pd.DataFrame(col)
    col_df = (col_df
              .rename(columns= {col_df.columns.values[0]: "Abundance"})
              .sort_values(by= ["Abundance"], ascending= False)
              .assign(cum_abundance= lambda x: x["Abundance"].cumsum()))
    col_df = col_df.loc[col_df['cum_abundance'] <= 95]
    proportion_col = col_df['Abundance'] / sum(col_df['Abundance'])
    shannon_col = proportion_col * np.log(proportion_col)
    shannon_index = -sum(shannon_col)
    richness = len(col_df)

    
    return(shannon_index, richness)

def apply_shannon95(fp):

    '''Takes in the filepath and calculates
    the shannon index from the dataframe given'''
    df = pd.read_csv(fp, sep= ',')
    shannon_df = (df
                  .iloc[:, 2:]
                  .apply(shannon_95, axis= 'index')
                  .transpose()
                  .rename(columns= {0: 'shannon_index',
                                    1: 'richness'}))
    return(shannon_df)

# run for genus healthy
genus_healthy95_path = Path('Outputs/healthy_df_out.csv')
genus_healthy95_shannon_df = apply_shannon95(genus_healthy95_path)
genus_healthy95_output_path = Path('graphs/diversity/genus_healthy95.csv')
genus_healthy95_shannon_df.to_csv(genus_healthy95_output_path)

# run for genus CRC
genus_CRC95_path = Path('Outputs/Stage_I_II_df_out.csv')
genus_CRC95_shannon_df = apply_shannon95(genus_CRC95_path)
genus_CRC95_output_path = Path('graphs/diversity/genus_CRC95.csv')
genus_CRC95_shannon_df.to_csv(genus_CRC95_output_path)







