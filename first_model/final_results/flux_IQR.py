####################
# flux_IQR.py
####################

'''Calculates the IQR (interquartile range), range
and median for all the flux values in variation_CRC
Example file is: variation_EX__4abut_medium.csv'''

import pandas as pd
import numpy as np
from pathlib import Path
import re

def summary_calc(fp, name):

    '''calculates the summary values range,
    IQR and median for the 10 samples in each
    dataframe'''

    
    met_df = pd.read_csv(fp, sep= ',')
    # don't want average as it will effect summary calcs
    met_df = met_df.loc[met_df['sample'] != "average"]
    met_df['flux'] = met_df['flux'] * -1
    maximum = met_df['flux'].max()
    minimum = met_df['flux'].min()
    median = met_df['flux'].median()
    Q1 = met_df['flux'].quantile(0.25)
    Q3 = met_df['flux'].quantile(0.75)
    IQR = Q3 - Q1
    # converts the summary data to a dataframe
    summary_data = {'max': maximum, 
                    'min': minimum, 
                    'median': median, 
                    'IQR': IQR}
    summary_df = pd.DataFrame(data = summary_data, index= [name])
    return summary_df

def summary_formatter(directory):

    '''Calculates the summary information
    for each file in the directory and
    puts into on summary dataframe'''
    
    count = 0
    for file in directory.glob('*csv'):
        crc_stem = file.stem
        if crc_stem.startswith('variation'):
            met_name = re.split(r'__|_', crc_stem)[2]
            # creates the summary dataframes and makes sure
            # there all concatanated together
            # itertaively adds each row to all summary df
            if count == 0:
                all_summary_df = summary_calc(file, met_name)
                count = 1
            else:
                summ_df = summary_calc(file, met_name)
                all_summary_df = pd.concat([all_summary_df, summ_df], axis = 0)
    return all_summary_df

# runs the main code
CRC_directory = Path('variation_CRC')
crc_summary_fp = CRC_directory / 'crc_var_summary.csv'
crc_summary_df = summary_formatter(CRC_directory)
crc_summary_df.to_csv(crc_summary_fp)

healthy_directory = Path('variation_healthy')
healthy_summary_fp = healthy_directory / 'healthy_var_summary.csv'
healthy_summary_df = summary_formatter(healthy_directory)
healthy_summary_df.to_csv(healthy_summary_fp)
