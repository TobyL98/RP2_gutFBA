######################
# sensitivity.py
######################

'''Performs a sensitivity analysis of
how percentage of the genera with a 
cretain abundance effect the community
modelling and community biomass generation'''

import pandas as pd
from pathlib import Path 
from average_abundance import average_abundance
from analysis import abundance_dict
from set_up import get_models

def main():

    '''Code will run through different cut off values of the
    percentage of abundance we want the genera to cover.
    Calculates the correct abundance and then runs the 
    analysis code to perform FBA'''

    # required inputs
    all_models_fp = Path('models')
    run_models_fp = Path('models_to_run_2')
    diet_medium_fp = Path("diet_info/average_EU_fluxes.tsv")
    overall_output_fp = Path("Results/sensitivity")
    matlab = False


    # list of cut offs to run
    cut_off_list = 0.3, 0.4, 0.5

    # read in files
    healthy_path = Path("Outputs/healthy_df_out.csv")
    healthy_df = pd.read_csv(healthy_path, sep = ',')

    Stage_I_II_path = Path("Outputs/Stage_I_II_df_out.csv")
    Stage_I_II_df = pd.read_csv(Stage_I_II_path, sep = ',')

    # loop through different cut offs
    for cut_off in cut_off_list:
        healthy_average_df = average_abundance(healthy_df, cut_off)
        # StageI_II_average_df = average_abundance(Stage_I_II_df, cut_off)



    # run analysis
    # getting correct models
    get_models(healthy_average_df, all_models_fp, run_models_fp, matlab)


    # Converting abundances from dataframe to dictionary
    abund_dict = abundance_dict(healthy_average_df, run_models_fp, matlab)

    # creating the community model
    com_model_obj = model_creation(run_models_fp, matlab)

    # generate the correct diet
    if diet_medium is not None:
        com_model_obj = generate_medium(diet_medium_fp, com_model_obj, output_fold_fp)


if __name__ == "__main__":
    main()