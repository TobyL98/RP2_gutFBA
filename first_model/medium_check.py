####################
# medium_check.py
####################

'''check whether the medium lower bounds correctly
match the maximum uptake flux for diet conditions'''

import pandas as pd 
import re

def id_sort(reaction_id):
    reaction_id_parts = reaction_id.split("_")
    if reaction_id_parts[1] == "":
        del(reaction_id_parts)[1]
    del(reaction_id_parts)[-1]
    reaction_id = "_".join(reaction_id_parts)
    return (reaction_id)

medium_flux_df = pd.read_csv("Results/medium_fluxes.csv", sep = ",", names = ["Reaction", "full_reaction","lower_bound", "upper_bound"])
western_flux_df = pd.read_csv("diet_info/western_fluxes.tsv", sep = "\t")

# sorting medium_flux_df
medium_flux_df["lower_bound"] = medium_flux_df["lower_bound"].apply(lambda bound: bound.split("-")[1])
medium_flux_df = medium_flux_df.loc[medium_flux_df["lower_bound"] != "1e"]
medium_flux_df = medium_flux_df.astype({'lower_bound': 'float'})
medium_flux_df.loc[:, 'Reaction'] = medium_flux_df['Reaction'].apply(id_sort)

#sorting western_flux_df
western_flux_df.loc[:, 'Reaction'] = western_flux_df['Reaction'].apply(lambda reaction: re.split(r"\[|\(", reaction)[0])

merged_flux_df = pd.merge(medium_flux_df, western_flux_df, on = ['Reaction'], how = "outer")

merged_flux_df['Correct'] = merged_flux_df['Flux Value'] == merged_flux_df['lower_bound']

print(merged_flux_df.head())
merged_flux_df.to_csv("medium_check.csv")