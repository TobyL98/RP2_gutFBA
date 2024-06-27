####################
# Abundance calculator
####################

# function will calculate
# average abundance
# and then filter by x number of species

import pandas as pd
from pathlib import Path

filepath = Path("top_25_species/CRC.csv")
df = pd.read_csv(filepath)
df = df.set_index('Species')

df["average_abundance"] = df.mean(axis = 1)

df = df.sort_values(by = 'average_abundance', ascending= False)
top10_df = df.head(10)
top10_df = top10_df['average_abundance']
top10_df = top10_df / sum(top10_df)
print(top10_df)



