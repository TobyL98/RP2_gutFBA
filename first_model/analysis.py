##################
# analysis.py
##################

'''Code uses single organism models to generate a community metabolic model
using PyCoMO and then runs FBA on the model. Outputs the Uptake fluxes and 
secretion fluxes into the medium. Also runs FVA to calculate the metabolite exchanges'''

import pandas as pd
import cobra
from pathlib import Path
import re
import time
import argparse
import sys
from set_up import get_models
from scipy import io


# import PyCoMo
# path_root = Path("C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research Project 2/Pycomo/PyCoMo/src")
path_root = Path("/Home/Documents/Toby/PyCoMo/src/pycomo")
sys.path.append(str(path_root))
import pycomo

# tests if an input file existsython 
def file_test(arg):
    p = Path(arg)
    if p.is_file():
        return p
    else:
        raise argparse.ArgumentTypeError("Argument {0} should be a valid file".format(arg))

# test if an input directory exists
def directory_test(arg):
    p = Path(arg)
    if p.is_dir():
        return p
    else:
        raise argparse.ArgumentTypeError("The input directory does not exist {0}". format(p))
    

def parse_args():

    description = """This code runs fixed abundance community metabolic model analysis using PyCoMo.
     The outputs are csv files of the uptake fluxes,
     secretion fluxes and the metabolic exchanges"""
    # set up argparse
    parser = argparse.ArgumentParser(description= description)

    # arguments
    parser.add_argument("-am", "--all_models",
                        help = """The directory containing every possible SBML model that could be required for the community model. 
                        This will then be filtered to just the models in the community model using the abundance dataframe. 
                        This is not required if you already have a folder with just the SBML models you need to run the community modelling.""",
                        type = directory_test)
    parser.add_argument("-rm", "--run_models",
                        help = """Directory containing only the community SBML models to run in the PyCoMO community modelling. Only an empty directory
                        is required if you are filtering using the data from the all_models directory""",
                        default = 'models_to_run',
                        type = directory_test)
    parser.add_argument("-a", "--abundance",
                        help = """The abundance dataframe with the abundances
                        for each organism in the model. Make sure the bacteria
                        are the same as in the model""",
                        default= 'Outputs/average/healthy_df_out_ave.csv',
                        type = file_test)
    parser.add_argument("-o", "--output_fold",
                        help = """The folder the outputs are saved into. It will create
                        three files containing met_exchanges, secretion fluxes and Uptake fluxes.""",
                        default = 'Results',
                        type = directory_test)
    parser.add_argument("-me", "--metExchange",
                        help = """Option for calculating the metabolic exchanges by running
                        FVA.""",
                        default= False,
                        action = 'store_true')
    parser.add_argument("-mat", "--matlab",
                        help = """Option for using matlab file (.mat) instead of the SBML format""",
                        default = False,
                        action = 'store_true')
    args = parser.parse_args()
    return(args)


def abundance_dict(df_path, FBA_models_path, matlab):
    '''Creates the abundance dictionary from
    an abundance dataframe
    INPUT: abundance dataframe'''

    print("\n#############################")
    print("Step 1: Creating the abundance dictionary")
    print("#############################")

    d_type = {"average_abundance": float}
    df = pd.read_csv(df_path, sep = ",", dtype = d_type)
    abun_df = df.set_index('Genus')['average_abundance']
    abun_dict = abun_df.to_dict()

    # checking model available for all bacteria in dictionary
    # if not will remove from dictionary
    if matlab == True:
        FBA_model_files = FBA_models_path.glob('*.mat')
    else:
        FBA_model_files = FBA_models_path.glob('*.xml')

    # get bacteria names from model
    species_we_have = []
    for model in FBA_model_files:

        species_name = model.stem
        species_we_have.append(species_name)
    
    # remove any bacteria from dictionary
    # with models not available
    for name in list(abun_dict.keys()):
        if name not in species_we_have:
            abun_dict.pop(name)
    
    # ensure all dictionary values sum to 1 after removal
    total_abun = sum(abun_dict.values())
    for key, value in abun_dict.items():
        new_value = value / total_abun
        abun_dict[key] = new_value

    return abun_dict

def model_creation(model_member_dir, matlab):

    '''Function loads the individual member models that will be used to create
    the community model and converts them to PyCoMo format. 
    Creates the community model object from these.'''

    print("\n#############################")
    print("Step 2: Creating the Community model")
    print("#############################")

    # load in the model members
    if matlab == True:
        named_models = pycomo.load_named_models_from_dir(model_member_dir, format = "mat")
    else:
        named_models = pycomo.load_named_models_from_dir(model_member_dir, format = "sbml")

    # check if all have a biomass objective function then it is working, will flag warning if not
    for model in named_models.values():
        objective = str(model.objective)
        try: 
            re.search(r"[Bb]iomass" , objective)
        except:
            print("""WARNING: Objective function not set to Biomass.
                  Therefore, code may fail""")
        del objective    

    # creating the single organism models in PyCoMo format
    single_org_models = []
    print("Creating single organism models for:")
    for name, model in named_models.items():
        print(name)
        single_org_model = pycomo.SingleOrganismModel(model, name)
        single_org_models.append(single_org_model)

        del single_org_model
    del named_models

    # creating the community model
    community_name = "gut_model"
    com_model_object = pycomo.CommunityModel(single_org_models, community_name)
    del single_org_models

    return com_model_object


def fixed_abundance(com_model_object, abundance_dict, output_folder):

    '''Function generates the fixed abundance model.
    Using generated abundance data '''

    print("\n#############################")
    print("\nStep 3: Creating Fixed abundance community model\n\n")
    print("#############################")

    # creating the abundance_dict
    abundance_total = sum(abundance_dict.values())
    if round(abundance_total, 3) != 1.00:
        print("""WARNING: abundance dictionary don't add to 1.
              Model may not work""")
        
    # converting to fixed abundance model and then adding the abundances from abundance_dict
    com_model_object.convert_to_fixed_abundance()
    com_model_object.apply_fixed_abundance(abundance_dict)

    # outputting abundance dictionary for results
    output_path = output_folder / r"abundance_dict.txt"
    with open(output_path, 'w') as f:  
        for key, value in abundance_dict.items():  
            f.write("{0}:{1}\n".format(key, value))
    f.close()


    print("\n\nFinished Creating Community Model!!!")
    return(com_model_object)


def final_analysis(com_model_object, output_dir):

    '''Function analyses the fixed abundance models.
    Gets the objective value, checks if solution is infeasible,
    and outputs uptake and secretion fluxes for the medium'''
    
    print("\n#############################")
    print("Step 4: Creating Outputs")
    print("#############################")
    solution = com_model_object.model.optimize()
    summary = com_model_object.summary()

    # check if is infesiable
    status = solution.status
    if str(status) == "infeasible":
        print("Solution infeasible. will not output results")
    else: 

        up_flux = summary.uptake_flux # how to get the uptake flux
        sec_flux = summary.secretion_flux # how to get the secretion flux

        up_flux_path = output_dir / r"up_flux.csv"
        sec_flux_path = output_dir / r"sec_flux.csv"
        up_flux.to_csv(up_flux_path)
        sec_flux.to_csv(sec_flux_path)


        # can get the objective value
        print("\nFlux of the objective value:")
        print(solution.objective_value) # how to get the objective_value

        print("""\nOutputs are:
              Uptake fluxes: {0}
              Secretion fluxes: {1}""".format(up_flux_path, sec_flux_path))

        return com_model_object

def metabolite_exchange(com_model_object, output_dir):

    '''Function uses FVA to calculate metabolite exchanges
    and cross feeding interactions'''

    print("\n#############################")
    print("Step 5: Calculating Metabolic Exchanges:")
    print("#############################")
    MES = com_model_object.potential_metabolite_exchanges()
    MES_path = output_dir / r"mes_results.csv"
    MES.to_csv(MES_path)
    del MES
    print("\n OUTPUT: {0}".format(MES_path))

def main():

    """Function that runs the code"""

    args = parse_args()

    # getting the models required for the specific Community
    # only runs if all_models input given
    all_models = args.all_models
    if all_models is not None:
        get_models(args.abundance, args.all_models, args.run_models, args.matlab)

    # Converting abundances from dataframe to dictionary
    abund_dict = abundance_dict(args.abundance, args.run_models, args.matlab)

    start_time = time.time()

    # creating the community model
    com_model_obj = model_creation(args.run_models, args.matlab)

    # switching to fixed abundance model and running
    com_model_obj = fixed_abundance(com_model_obj, abund_dict, args.output_fold)

    # Doing the analysis to obtain outputs
    com_model_obj_final = final_analysis(com_model_obj, args.output_fold)

    # Working out possible metabolic exchnages
    metab_exchange = args.metExchange
    if metab_exchange == True and com_model_obj_final:
        metabolite_exchange(com_model_obj_final, args.output_fold)
    elif metab_exchange == False:
        print("\nNOT calculating metabolic exchanges")
    else:
        print("""\nNOTICE: Not calculating metabolite exchanges as
              Solution may be infeasible""")

    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60

    del com_model_obj

    print("\nTotal running time (minutes):{0}".format(elapsed_time))

if __name__ == "__main__":
    main()
