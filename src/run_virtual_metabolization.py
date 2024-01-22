import pandas as pd
import numpy as np
import math
from subprocess import call
import subprocess
from rdkit.Chem import SaltRemover, MolStandardize
from molvs import Standardizer
from rdkit import Chem
from rdkit import RDLogger
import shutil
import traceback
import requests
import pkg_resources
import zipfile

lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

import os
import sys

# Function to check and install packages
def install_package(package):
    try:
        pkg_resources.require(package)
        print(f"{package} is already installed.")
    except pkg_resources.DistributionNotFound:
        print(f"{package} not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Import and install SyGMa
install_package("sygma")

# Import SyGMa after ensuring it's installed
import sygma

#Write out the results
def export_for_SIRIUS(input_df_name, compound_name):
    df = pd.read_csv(input_df_name, sep='\t')
    print('==== Export for SIRIUS ====')
    print('Number of metabolites = ' +str(df.shape[0]))
       
    try:
        #SyGMa
        df_new = df.dropna(subset=['metabolite']) 
        df_new = df_new[['metabolite', compound_name+'_SyGMa']]
        df_new = df_new.rename(columns={"metabolite": "Smiles", compound_name+'_SyGMa': "name"})
    except:
        #BioTransformer
        df = df.dropna(subset=['SMILES']) 
        df_new = df.copy()
        df_new.loc[:, 'Biotransformed ' + compound_name] = df_new['Reaction'].astype(str).apply(lambda x: x.split(';')[0]) + " > " + df_new[compound_name].astype(str)
        df_new = df_new[['SMILES', 'Biotransformed '+compound_name]]
        df_new = df_new.rename(columns={"SMILES": "Smiles", 'Biotransformed '+compound_name: "name"})

    df_new = df_new.drop_duplicates(subset='name', keep='first')
    df_new = df_new.drop_duplicates(subset='Smiles', keep='first')

    print('Number of unique metabolites considered = ' +str(df_new.shape[0])) 
    df_new.to_csv(input_df_name[:-4]+'_SIRIUS.tsv', sep = '\t', index = False)
    

def export_for_NAP(input_df_name, compound_name):
    df = pd.read_csv(input_df_name, sep='\t')
    print('==== Export for NAP ====')
    print('Number of metabolites = ' +str(df.shape[0]))  

    try:
        #SyGMa
        df_new = df.dropna(subset=['metabolite']) 
        df_new = df_new[['metabolite', compound_name+'_SyGMa']]
        df_new = df_new.drop_duplicates(subset='metabolite', keep='first')
        df_new = df_new.drop_duplicates(subset=compound_name+'_SyGMa', keep='first')
    except:
        #BioTransformer
        df_new = df.copy()
        df_new = df_new.dropna(subset=['SMILES'])
        df_new.loc[:, 'Biotransformed ' + compound_name] = df_new['Reaction'].astype(str).apply(lambda x: x.split(';')[0]) + " > " + df_new[compound_name].astype(str)
        df_new = df_new[['SMILES', 'Biotransformed '+compound_name]]
        df_new = df_new.drop_duplicates(subset='SMILES', keep='first')
        #df_new = df_new.drop_duplicates(subset='SMILES',compound_name+'_BioTransformer', keep='first')


    print('Number of unique metabolites considered = ' +str(df_new.shape[0]))   
    df_new.to_csv(input_df_name[:-4]+'_NAP.tsv', sep = '\t', index = False, header= False)


# Run SyGMa on two lists of SMILES and compound name in batch
def run_sygma_batch(list_smiles, list_compound_name, phase_1_cycle, phase_2_cycle, top_sygma_candidates, output_name, compound_name):
    
    #Define the batch_size based on number of phase
    batch_size = int(100/(math.exp(phase_1_cycle))/(math.exp(phase_2_cycle)))
    
    print('#### Starting SyGMa computation')
    print('Number of compounds that will be subject to SyGMa metabolisation = '+str(len(list_smiles)))
    print('Number of phase I cycles = '+str(phase_1_cycle,)) 
    print('Number of phase II cycles = '+str(phase_2_cycle))
    print('Keep only the top SyGMa candidates = '+str(top_sygma_candidates))
    print('If you are running many compounds or cycles, and maxing out RAM memory available, you can decrease the batch size. Otherwise the value can be increased for faster computation.')
    print('======')
    print('Please wait')
    print('======')

    # Define SyGMa parameters      
    scenario = sygma.Scenario([
        [sygma.ruleset['phase1'], int(phase_1_cycle)],
        [sygma.ruleset['phase2'], int(phase_2_cycle)]
    ])
    
    # Create placeholder dataframe file
    df_master = pd.DataFrame(columns=['metabolite', 'score', 'parent', 'pathway', 'Compound_Name'])

    # Creating batches
    for i in range(0, len(list_smiles), batch_size):
        list_smiles_batch = list_smiles[i:i + batch_size]
        list_compound_name_batch = list_compound_name[i:i + batch_size]

        # Running SyGMa on each batch
        for a, b in zip(list_smiles_batch, list_compound_name_batch):
            try:
                mol = Chem.MolFromSmiles(a)
                metabolic_tree = scenario.run(mol)
                metabolic_tree.calc_scores()

                metabolites = metabolic_tree.to_smiles()
                metabolite_info = metabolic_tree.to_list()

                # Extract metabolite data and create a DataFrame
                metabolite_data = [{
                    'metabolite': met[0],
                    'score': round(met[1], 3) if met[1] != 'N/A' else 'N/A',
                    'parent': metabolites[0][0],
                    'pathway': info['SyGMa_pathway'][:200],
                    'Compound_Name': b
                } for met, info in zip(metabolites, metabolite_info)]

                # Limit to top candidates and remove all-NaN columns
                if metabolite_data:
                    df = pd.DataFrame(metabolite_data)
                    df = df.head(top_sygma_candidates + 1)  # +1 to include the parent metabolite
                    df = df.dropna(axis=1, how='all')
                    df_master = pd.concat([df_master, df], ignore_index=True)

            except Exception as e:
                print(f"Error processing molecule: {a}. Error: {e}")
                traceback.print_exc()  # Print the full traceback

        print(f'Batch {i // batch_size + 1}/{math.ceil(len(list_smiles) / batch_size)} completed')

    df_master['score'] = df_master['score'].round(3)
    df_master[compound_name+'_SyGMa'] = df_master['pathway'] + df_master['score'].astype(str) + '; ' + df_master[compound_name]
    df_master[compound_name+'_SyGMa'] = df_master[compound_name+'_SyGMa'].str.replace("\n", "")

    print(f'Number of SyGMA candidates = {df_master.shape[0]}\nNumber of unique SyGMA candidates = {len(df_master.metabolite.unique())}\n=== SyGMa COMPLETED ===')
    
    df_master.to_csv(str(output_name)+'_results_vm_NAP_SyGMa.tsv', sep='\t', index=False)
    run_sygma_batch.df_master = df_master

    file_name = f"{output_name}_results_vm_NAP_SyGMa"
    run_sygma_batch.file_name_sygma = f"{file_name}.tsv"

    print('  ')
    export_for_NAP(run_sygma_batch.file_name_sygma, 'Compound_Name')
    print('  ')
    export_for_SIRIUS(run_sygma_batch.file_name_sygma, 'Compound_Name')
    print('  ')
    print('---- DOWNLOAD THE RESULTS BELOW ----')
    run_sygma_batch.file_name_sygma_nap = f"{file_name}_NAP.tsv"
    run_sygma_batch.file_name_sygma_sirius = f"{file_name}_SIRIUS.tsv"
    run_sygma_batch.markdown_link_sygma = f"[View/Download the vm_NAP SyGMa results from {output_name}](./{run_sygma_batch.file_name_sygma})."
    run_sygma_batch.markdown_link_sygma_nap = f"[View/Download the vm_NAP SyGMa results for Network Annotation Propagation from {output_name}](./{run_sygma_batch.file_name_sygma_nap})."
    run_sygma_batch.markdown_link_sygma_sirius = f"[View/Download the vm_NAP SyGMa results for SIRIUS from {output_name}](./{run_sygma_batch.file_name_sygma_sirius})."

#Defining a function to check the Java version
def check_java_version():
    """Check and print the installed Java version."""
    try:
        java_version = subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT)
        print(java_version.decode())
        print("Java version check completed.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while checking Java version: {e.output.decode()}")
    except Exception as e:
        print(f"An unexpected error occurred while checking Java version: {e}")

#Defining a function to download and unzip Biotransformer
def old_download_and_unzip_biotransformer():
    """Download and unzip the BioTransformer software if not already present."""

    if not os.path.exists('biotransformer3.zip'):
        print("Downloading biotransformer3 ...")
        url = 'https://bitbucket.org/wishartlab/biotransformer3.0jar/get/6432cf887ed70.zip'
        try:
            response = requests.get(url, allow_redirects=True)
            response.raise_for_status()
            with open('biotransformer3.zip', 'wb') as file:
                file.write(response.content)
        except requests.RequestException as e:
            print(f"An error occurred while downloading BioTransformer: {e}")
            return
        print("BioTransformer downloaded.")
    else:
        print("BioTransformer was already downloaded - skipping download.")

    if not os.path.exists('BioTransformer3.0_20230525.jar'):
        try:
            with zipfile.ZipFile('biotransformer3.zip', 'r') as zip_ref:
                zip_ref.extractall('.')
            print("BioTransformer is unzipped.")
        except Exception as e:
            print(f"An error occurred while unzipping BioTransformer: {e}")
    else:
        print("BioTransformer is already unzipped - skipping unzip.")

def download_and_unzip_biotransformer():
    """Download and unzip the BioTransformer software if not already present."""

    url = 'https://bitbucket.org/wishartlab/biotransformer3.0jar/get/6432cf887ed70.zip'
    local_zip_file = 'biotransformer3.zip'
    local_jar_file = 'BioTransformer3.0_20230525.jar'

    if not os.path.exists(local_zip_file):
        print(f"Downloading {local_zip_file} from {url}...")
        try:
            response = requests.get(url, allow_redirects=True)
            if response.status_code == 200:
                with open(local_zip_file, 'wb') as file:
                    file.write(response.content)
                print(f"Downloaded {local_zip_file}.")
            else:
                print(f"Failed to download. Status code: {response.status_code}")
        except requests.RequestException as e:
            print(f"An error occurred while downloading BioTransformer: {e}")
            return
    else:
        print(f"{local_zip_file} was already downloaded - skipping download.")

    if not os.path.exists(local_jar_file):
        print(f"Unzipping {local_zip_file}.")
        try:
            with zipfile.ZipFile(local_zip_file, 'r') as zip_ref:
                zip_ref.extractall('.')
            print(f"{local_zip_file} is unzipped.")
        except Exception as e:
            print(f"An error occurred while unzipping BioTransformer: {e}")
    else:
        print(f"{local_jar_file} is already unzipped - skipping unzip.")

#Prepare the execution environment for biotransformer3 by creating or cleaning the output folder
def prepare_environment(source_dir, dest_dir):
    """Prepare the working environment for biotransformer3 by moving files and cleaning up directories."""
    if os.path.exists(source_dir):
        readme_path = os.path.join(source_dir, "README.md")
        new_readme_path = os.path.join(source_dir, "README_BioTransformer.md")
        if os.path.exists(readme_path):
            os.rename(readme_path, new_readme_path)
            print("Renamed README.md to README_BioTransformer.md")
            for filename in os.listdir(source_dir):
                source_file = os.path.join(source_dir, filename)
                dest_file = os.path.join(dest_dir, filename)
                shutil.move(source_file, dest_file)
    else:
        print(f"Directory {source_dir} does not exist.")

#Defining a function to check if the biotransformation type is valid
def validate_biotransformation_type(type_of_biotransformation):
    list_of_biotransformation = ['ecbased', 'cyp450', 'phaseII', 'hgut', 'superbio', 'allHuman', 'envimicro']
    if type_of_biotransformation in list_of_biotransformation:
        print('     Biotransformation: '+type_of_biotransformation)
        return True
    else:
        print('Check the type/spelling of the biotransformation!')
        return False

#Defining a function to create or clear the output folder        
def create_or_clear_output_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    else:
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f'Failed to delete {file_path}. Reason: {e}')

def prepare_for_bio3(type_of_biotransformation, list_smiles):
    """Prepare the working environment for biotransformer3 by checking Java, Downloading and unziping BioTransformer3, prepare the folder, validate the biotransformation type, creating necessary files and cleaning up directories."""

    output_folder = "biotransformer_results"

    # Check Java version
    try:
        check_java_version()
        print("Java version check completed.")
    except Exception as e:
        print(f"Error checking Java version: {e}")

    # Download and unzip BioTransformer
    try:
        download_and_unzip_biotransformer()
        print("BioTransformer download and unzip completed.")
    except Exception as e:
        print(f"Error downloading and unzipping BioTransformer: {e}")

    # Prepare the working environment
    try:
        prepare_environment('wishartlab-biotransformer3.0jar-6432cf887ed7', '.',)
        print("BioTransformer environment preparation completed.")
    except Exception as e:
        print(f"Error preparing the working environment: {e}")

    # Validate the biotransformation type
    try:
        validate_biotransformation_type(type_of_biotransformation)
        print("Biotransformation type validation completed.")
    except Exception as e:
        print(f"Error validating biotransformation type: {e}")

    # Create the output folder if it doesn't exist, clear if it does
    try:
        create_or_clear_output_folder(output_folder)
        print("Output folder creation/clearing completed.")
    except Exception as e:
        print(f"Error creating or clearing the output folder: {e}")

    try:    
        os.remove("biotransformer_temp_output.csv")
    except FileNotFoundError:
        pass

    print('######  Parameters checking for Biotransformer3 completed. Now starting the computation...')
    print('######  Running BioTransformer takes approximatively 3-60 secs per compound depending on biotransformation parameters')
    print('     Number of compounds being virtually metabolized with BioTransformer =  '+str(len(list_smiles)))

# Run BioTransformer3 on two lists of SMILES and compound name - 2211
def run_biotransformer3(mode, list_smiles, list_compound_name, type_of_biotransformation, number_of_steps, output_name):
    """Run BioTransformer3 on SMILES string"""

    output_folder = "biotransformer_results"
    col_list = ['InChIKey',	'SMILES','PUBCHEM_CID','Molecular formula','Major Isotope Mass',
                    'Metabolite ID','cdk:Title','Reaction','Reaction ID', 'Enzyme(s)','Biosystem','Precursor SMILES',
                    'Precursor Major Isotope Mass']
    
    if not os.path.exists('BioTransformer3.0_20230525.jar'):
        print("Error: BioTransformer3.0_20230525.jar not found. Please run prepare_for_bio3 first.")
        return

    counter = 0
    for a, b in zip(list_smiles, list_compound_name):
            
        counter +=1
        output_filename = f"{output_folder}/biotransformer_output_{counter}.csv"
            
        try:
            if mode == 'standard' or mode == 'btType':
                biotransformcall = 'java -jar BioTransformer3.0_20230525.jar -k pred -b ' + type_of_biotransformation +' -ismi "' + a +'" -ocsv '+str(output_filename)+' -s '+str(number_of_steps)
            elif mode.str.startswith('-k'):
                biotransformcall = 'java -jar BioTransformer3.0_20230525.jar -ismi "' + a +'" -ocsv '+str(output_filename)+' '+mode

            biotransformcall = biotransformcall.split() # because call takes a list of strings 
            call(biotransformcall)

            try:
                # Check if the output file exists and is not empty
                if os.path.exists(output_filename) and os.path.getsize(output_filename) > 0:
                    df_bio = pd.read_csv(output_filename, sep=',', usecols=col_list)
                    df_bio['Parent_Compound_Name'] = b
                    #df_bio['Biotransformed_Compound_Name'] = str(df_bio['Reaction'][:25]) + str(df_bio['Parent_Compound_Name'])
                    df_bio.to_csv(output_filename, sep='\t', index=False)
                else:
                    print(f"No output for compound {b}")

            except Exception as e:
                print(f"Error with BioTransformer for compound {b}: {e}")

        except:
            print('          ! No candidate or error with BioTransformer for compound n'+str(counter)+' - it will be ignored')
            print('          '+b)
            print('          For:  '+a+' . Check the SMILES on http://biotransformer.ca and/or the GNPS library entry.')
            pass
    
    dataframes = []
    for i in range(1, counter + 1):
        temp_output_filename = f"{output_folder}/biotransformer_output_{i}.csv"
        if os.path.exists(temp_output_filename):
            temp_df = pd.read_csv(temp_output_filename, sep='\t')
            dataframes.append(temp_df)

    # Concatenate all DataFrames in the list
    df2_bio = pd.concat(dataframes, ignore_index=True)

    print(' - ')
    print('####### VIRTUAL METABOLISATION COMPLETED FOR ALL COMPOUNDS #######')
    print('Total number of BioTransformer candidates = '+str(df2_bio.shape[0])+' - and '+str(len(df2_bio.SMILES.unique()))+" are unique")
    #print('Total number of unique BioTransformer candidates = '+str(len(df2_bio.SMILES.unique())))

    file_name = f"{output_name}_results_vm_NAP_BioTransformer"
    run_biotransformer3.file_name_biotransf = f"{file_name}.tsv"
    run_biotransformer3.file_name_biotransf_nap = f"{file_name}_NAP.tsv"
    run_biotransformer3.file_name_biotransf_sirius = f"{file_name}_SIRIUS.tsv"

    df2_bio.to_csv(run_biotransformer3.file_name_biotransf, sep='\t', index = False)
    
    print('  ')
    print('  ')
    export_for_NAP(run_biotransformer3.file_name_biotransf, 'Parent_Compound_Name')
    print('  ')
    export_for_SIRIUS(run_biotransformer3.file_name_biotransf, 'Parent_Compound_Name')
    print('  ')
    print('##### DOWNLOAD THE RESULTS BELOW #####')
    run_biotransformer3.markdown_link_biotransf = f"[View/Download the vm_NAP results from {output_name}](./{run_biotransformer3.file_name_biotransf})."
    run_biotransformer3.markdown_link_biotransf_nap = f"[View/Download the vm_NAP BioTransformer results for Network Annotation Propagation from {output_name}](./{run_biotransformer3.file_name_biotransf_nap})."
    run_biotransformer3.markdown_link_biotransf_sirius = f"[View/Download the vm_NAP BioTransformer results for SIRIUS from {output_name}](./{run_biotransformer3.file_name_biotransf_sirius})."