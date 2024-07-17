######################
# sensitivity.py
######################

'''Performs a sensitivity analysis of
how percentage of the genera with a 
cretain abundance effect the community
modelling and community biomass generation'''

import pandas as pd
from pathlib import Path 
import time
from average_abundance import average_abundance
from analysis import abundance_dict, model_creation, generate_medium, fixed_abundance, final_analysis
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
    overall_output_fp = Path("Results/sensitivity_CRC")
    matlab = False


    # list of cut offs to run
    cut_off_list = [0.95]

    # read in files
    healthy_path = Path("Outputs/healthy_df_out.csv")
    healthy_df = pd.read_csv(healthy_path, sep = ',')

    Stage_I_II_path = Path("Outputs/Stage_I_II_df_out.csv")
    Stage_I_II_df = pd.read_csv(Stage_I_II_path, sep = ',')

    # loop through different cut offs
    objective_flux_list = []
    for cut_off in cut_off_list:
        print("\n###############")
        print("Running FBA of abundance cut off {0}".format(cut_off))
        print("###############")

        start_time = time.time()

        average_df = average_abundance(Stage_I_II_df, cut_off)

        # output folder
        cut_off_name = str(cut_off).split('.')[1]
        output_dir_path = overall_output_fp / r"cutoff_{0}".format(cut_off_name)
        print(output_dir_path)
        output_dir_path.mkdir()

        # run analysis
        # getting correct models
        get_models(average_df, all_models_fp, run_models_fp, matlab)


        # Converting abundances from dataframe to dictionary
        abund_dict = abundance_dict(average_df, run_models_fp, matlab)

        # creating the community model
        com_model_obj = model_creation(run_models_fp, matlab)

        # generate the correct diet
        if diet_medium_fp is not None:
            com_model_obj = generate_medium(diet_medium_fp, com_model_obj, output_dir_path)
            

        # creating fixed abundance model
        com_model_obj = fixed_abundance(com_model_obj, abund_dict, output_dir_path)

        # final analysis
        com_model_obj_final = final_analysis(com_model_obj, output_dir_path)

        end_time = time.time()
        elapsed_time = (end_time - start_time) / 60
        print("\nElapsed time for run: {0}".format(elapsed_time))
        time_filepath = output_dir_path / r"time_taken.txt"
        with open(time_filepath, "w") as f:
            f.write("Time taken to run:{0}".format(elapsed_time))
        f.close()
            
           
        

if __name__ == "__main__":
    main()
