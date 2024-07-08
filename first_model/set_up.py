##################
# set_up.py
##################

'''This code sets up the correct models ready for
analysis using PyCoMo in analysis.py'''

from pathlib import Path
import pandas as pd
import shutil
import re


def get_models(df_path, model_path, FBA_models_path, matlab):

    '''Function will get the models that match the
    genus in the inputted Dataframe
    Requirements: must have dataframe column 'Genus'''

    print("\n#############################")
    print("Pre-step: Obtaining the Correct Models")
    print("#############################")
    # deleting any exisiting models in FBA models
    # folder
    old_model_files = list(FBA_models_path.glob('*.xml')) + list(FBA_models_path.glob('*.mat'))
    for filepath in old_model_files:
        filepath.unlink()


    if isinstance(df_path, pd.DataFrame):
        df = df_path
    else:
        df = pd.read_csv(df_path, sep = ",")

    # SBML or matlab files of all metabolic models
    if matlab == True:
        model_files = model_path.glob('*.mat')
    else:
        model_files = model_path.glob('*.xml')
    

    for filepath in model_files:
        model_name = filepath.stem
        
        genus_name = str(model_name).split("_")[0]
        
        # if model name in the dataframe
        genus_df = df.loc[df['Genus'] == genus_name]
        if genus_df.shape[0] == 1:

            # will copy model to the destination folder
            # where FBA will be run
            if matlab == True:
                dest_path = FBA_models_path / "{0}.mat".format(genus_name) 
            else:
                dest_path = FBA_models_path / "{0}.xml".format(genus_name) 
            try:
                shutil.copy(filepath, dest_path)
            except shutil.SameFileError:
                print("Source and destination represents the same file.")
            except PermissionError:
                print("Permission denied.")
            except:
                print("Error occured while copying file")


def main():

    df_path = Path("Outputs/average/healthy_df_out_ave.csv")
    model_path = Path("Agora_Western/sbml")
    models_to_run_path = Path("models_to_run_2")
    matlab = False

    get_models(df_path, model_path, models_to_run_path, matlab)


if __name__ == "__main__":
    main()