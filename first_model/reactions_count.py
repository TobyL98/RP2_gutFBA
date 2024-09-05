###################
# reactions_model
###################

'''Counts the number of reactions for
all the models in the model_to_run folder'''

import cobra
from pathlib import Path
import math

folder_path = Path('models_to_run')
reactions_count_list = []
for model_path in folder_path.glob('*.xml'):
    model = cobra.io.read_sbml_model(model_path)
    reactions_count = len(model.reactions)
    model_name = model_path.stem
    print("{0}, Reaction count = {1}".format(model_name, reactions_count))
    reactions_count_list.append(reactions_count)

print(sum(reactions_count_list) / len(reactions_count_list))
