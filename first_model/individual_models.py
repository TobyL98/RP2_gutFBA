#####################
# individual_models.py
#####################

'''Runs the community model for each individual
member under average_european constraints'''

from pathlib import Path
import pandas as pd
import shutil 

def individual_models(models_dir,models_run_dir, diet_fp):
    '''Runs the community model for each individual
    member under average_european constraints'''

    for model_path in models_dir.iterdir():
        
        # deleting existing models
        for old_model_path in models_run_dir.iterdir():
            old_model_path.unlink()

        model_name = model_path.name
        dest_path = models_run_dir / model_name
        
        try:
            shutil.copy(model_path, dest_path)
        except shutil.SameFileError:
            print("Source and destination represents the same file.")
        except PermissionError:
            print("Permission denied.")
        except:
            print("Error occured while copying file")





models_path = Path('models_to_run')
models_run_path = Path('models_to_run_2')
diet_path = Path('diet_info/average_EU_fluxes.tsv')

individual_models(models_path, models_run_path, diet_path)
