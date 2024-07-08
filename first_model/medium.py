#######################
# medium.py
#######################

import pandas as pd
from pathlib import Path

def medium(diet_filepath, community_model_obj):
    '''
    Generates a community model specific medium
    from the reaction fluxes specified in a tsv
    file from the diet calculator in vmh.life
    https://www.vmh.life/#nutrition

    :param diet_filepath: the filepath for tsv file with the diet fluxes
    :param com_model: the community model object from PyCoMO
    '''
    parent_path = Path().absolute().parent
    diet_filepath = parent_path / r"diet_info/western_fluxes.tsv"
    diet_df = pd.read_csv(diet_filepath, sep = "\t")
    diet_df.loc[:, 'Reaction'] = diet_df['Reaction'].apply(lambda reaction: reaction
                                                           .replace("[", "(")
                                                           .replace("]", ")")
                                                           )
    diet_reactionID_list = list(diet_df.loc[:, 'Reaction'])

    # getting all the models in the community model
    model_members = community_model_obj.member_models

    #creating the medium with the diet values
    medium_exchange_dict = {}
    for single_org_model in model_members:
        model = single_org_model.model

        for reaction in model.exchanges:
            reaction_id = reaction.id
            if reaction_id in medium_exchange_dict.values():
                continue
            elif reaction_id in diet_reactionID_list:
                flux = diet_df.loc[diet_df['Reaction'] == reaction_id, 'Flux Value'].values
                medium_exchange_dict[reaction_id] = flux
            else:
                medium_exchange_dict[reaction_id] = 1e-6
    print(medium_exchange_dict["EX_strch1(e)"])