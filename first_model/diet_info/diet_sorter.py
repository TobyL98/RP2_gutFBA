######################
# diet_sorter.py
######################

import pandas as pd
from pathlib import Path

file_path = Path("fluxes.tsv")
diet_flux_df = pd.read_csv(file_path, sep = "\t")

print(diet_flux_df.head())


'''Already in the format I need so add thsi code elsewhere'''