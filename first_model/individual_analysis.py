import pandas as pd
import numpy as np
from pathlib import Path
import cobra 
import time
import sys
import random
from average_abundance_script import average_abundance
from analysis import abundance_dict, model_creation, generate_medium, fixed_abundance, final_analysis
from set_up import get_models

def extract_individual(df, FBA_models_path, cut_off, individual_num):
    
    '''Just gets an individuals bacterial abundances
    from a df with abundance results for many individuals.
    Then creates the dictionary which cane be usedf for analysis'''
    
    sample_df = df.iloc[: , [0, 1, individual_num]]

    # rename final column to abundance
    sample_df = sample_df.rename(columns = {sample_df.columns[2]: "abundance"})
    sample_df['abundance'] = sample_df['abundance'] / sum(sample_df['abundance'])

    # cumalative abundance
    sample_df = sample_df.sort_values(by = ['abundance'], ascending = False)
    sample_df['cumulative_abundance'] = sample_df['abundance'].cumsum()

    cut_off_df = sample_df.loc[sample_df['cumulative_abundance'] <= cut_off]

    return cut_off_df



def main():
    
    '''Runs the main code for analysis'''

    # required inputs
    input_path = Path("Outputs/healthy_df_out.csv")
    input_df = pd.read_csv(input_path, sep = ',')

    all_models_fp = Path("models")
    run_models_fp = Path("models_to_run_2")
    diet_medium_fp = Path("diet_info/average_EU_fluxes.tsv")
    overall_output_fp = Path("final_results/variation_healthy")
    matlab = False
    threshold = 0.95
    num_samples = 5

    # generating the random individuals to use
    np.random.seed(42)
    overall_samples_num = input_df.shape[1]
    # starts at 3 to avoid species and genus columns
    random_samples = np.random.randint(3, overall_samples_num, size=5)

    for sample_num in random_samples:
        start_time = time.time()
        print("Running sample{0}".format(sample_num))
        # get df of individuals abundance
        abundance_df = extract_individual(input_df, all_models_fp, threshold, sample_num)
        # output directory
        output_dir_path = overall_output_fp / r"sample_{0}".format(sample_num)

        # run analysis
        # getting correct models
        get_models(abundance_df, all_models_fp, run_models_fp, matlab)
        # Converting abundances from dataframe to dictionary
        abund_dict = abundance_dict(abundance_df, run_models_fp, matlab)
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