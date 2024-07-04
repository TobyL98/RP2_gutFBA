###################
# copy_models.py
###################

'''Deletes SBML models that don't match the model format
already created. This model format is the models in the
"models" folder'''

from pathlib import Path

def copy_models(correct_path, all_path, xml_mat):
    '''Copies SBML models from one folder to
another. Only copies models that match the
models in another folder'''

    correct_files = correct_path.glob('*.xml')

    if xml_mat == "mat":
        all_files = all_path.glob('*.mat')
    else:
        all_files = all_path.glob('*.xml')

    correct_filestem_list = []
    for correct_filepath in correct_files:
        correct_filestem = correct_filepath.stem
        correct_filestem_list.append(correct_filestem)
    
    for all_filepath in all_files:
        all_filestem = all_filepath.stem
       
        if all_filestem not in correct_filestem_list:
           all_filepath.unlink()





correct_models_path = Path('models')
all_models_path = Path('Agora_Western/mat')
xml_or_mat = "mat"

copy_models(correct_models_path, all_models_path, xml_or_mat)