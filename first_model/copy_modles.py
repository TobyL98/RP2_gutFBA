###################
# copy_models.py
###################

'''Deletes SBML models that don't match the model format
already created. This model format is the models in the
"models" folder'''

from pathlib import Path

def copy_models(correct_path, all_path):
    '''Copies SBML models from one folder to
another. Only copies models that match the
models in another folder'''

    correct_files = correct_path.glob('*.xml')
    all_files = all_path.glob('*.xml')

    correct_filename_list = []
    for correct_filepath in correct_files:
        correct_filename = correct_filepath.name
        correct_filename_list.append(correct_filename)
    
    for all_filepath in all_files:
        all_filename = all_filepath.name
       
        if all_filename not in correct_filename_list:
           all_filepath.unlink()





correct_models_path = Path('models')
all_models_path = Path('Agora_HighFiber/sbml')

copy_models(correct_models_path, all_models_path)