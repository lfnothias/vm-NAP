import pandas as pd
import numpy as np
import math
from subprocess import call
import subprocess
import sygma
from rdkit.Chem import SaltRemover, MolStandardize
from molvs import Standardizer
from rdkit import Chem
from rdkit import RDLogger
import shutil
import traceback

lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

import os
import sys

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
    
    df_master.to_csv(str(output_name)+'_results_vm-NAP_SyGMa.tsv', sep='\t', index=False)
    run_sygma_batch.df_master = df_master

    file_name = f"{output_name}_results_vm-NAP_SyGMa"
    run_sygma_batch.file_name_sygma = f"{file_name}.tsv"

    print('  ')
    export_for_NAP(run_sygma_batch.file_name_sygma, 'Compound_Name')
    print('  ')
    export_for_SIRIUS(run_sygma_batch.file_name_sygma, 'Compound_Name')
    print('  ')
    print('---- DOWNLOAD THE RESULTS BELOW ----')
    run_sygma_batch.file_name_sygma_nap = f"{file_name}_NAP.tsv"
    run_sygma_batch.file_name_sygma_sirius = f"{file_name}_SIRIUS.tsv"
    run_sygma_batch.markdown_link_sygma = f"[View/Download the vm-NAP SyGMa results from {output_name}](./{run_sygma_batch.file_name_sygma})."
    run_sygma_batch.markdown_link_sygma_nap = f"[View/Download the vm-NAP SyGMa results for Network Annotation Propagation from {output_name}](./{run_sygma_batch.file_name_sygma_nap})."
    run_sygma_batch.markdown_link_sygma_sirius = f"[View/Download the vm-NAP SyGMa results for SIRIUS from {output_name}](./{run_sygma_batch.file_name_sygma_sirius})."


# Run BioTransformer3 on two lists of SMILES and compound name - 2211
def run_biotransformer3(mode, list_smiles, list_compound_name, type_of_biotransformation, number_of_steps, output_name):
    
    try:
        # Check Java version
        java_version = subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT)
        print(java_version.decode())

        # Check if the JAR file exists, and download it if it doesn't
        if not os.path.exists('biotransformer3.zip'):
            print(f"Downloading biotransformer3 ...")
            subprocess.check_call(["curl", "-C", "-", "https://bitbucket.org/wishartlab/biotransformer3.0jar/get/6432cf887ed70.zip", "-o", 'biotransformer3.zip'])
            print("Download complete.")
        else:
            print(f"BioTransformer was already downloaded - skipping download.")
        
        try:
            if not os.path.exists('BioTransformer3.0_20230525.jar'):
                subprocess.check_call(["unzip", "-o", "-q", "biotransformer3.zip"])
            else:
                print(f"BioTransformer is already unzipped - skipping unzip.")

            print(" ")
            source_dir = 'wishartlab-biotransformer3.0jar-6432cf887ed7'
            dest_dir = '.'
            # Check if the source directory exists
            if os.path.exists(source_dir):
                # Move each file in the source directory to the destination directory
                for filename in os.listdir(source_dir):
                    source_file = os.path.join(source_dir, filename)
                    dest_file = os.path.join(dest_dir, filename)
                    shutil.move(source_file, dest_file)
            else:
                print(f"Directory {source_dir} does not exist.")

        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            
        # Prepare pandas tables
        df_bio = pd.DataFrame(data=None)
        print('######  Running BioTransformer takes approximatively 3-60 secs per compounds depending on biotransformation parameters')
        print('     Number of compounds being virtually metabolized with BioTransformer =  '+str(len(list_smiles)))
        
        col_list = ['InChIKey',	'SMILES','PUBCHEM_CID','Molecular formula','Major Isotope Mass',
                    'Metabolite ID','cdk:Title','Reaction','Reaction ID', 'Enzyme(s)','Biosystem','Precursor SMILES',
                    'Precursor Major Isotope Mass']


        list_of_biotransformation = ['ecbased','cyp450','phaseII','hgut','superbio','allHuman','envimicro']
        
        if type_of_biotransformation not in list_of_biotransformation:
            print('Check the type/spelling of the biotransformation !!!')
            return None
        elif type_of_biotransformation in list_of_biotransformation:
            print('     Biotransformation: '+type_of_biotransformation)
            print('     Please wait for the computation ...')


        output_folder = "biotransformer_results"

        # Create the output folder if it doesn't exist
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        else:
            # If the folder exists, clear its contents
            for filename in os.listdir(output_folder):
                file_path = os.path.join(output_folder, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print(f'Failed to delete {file_path}. Reason: {e}')

        counter = 0
            
        # Iterative into the lists and run BioTransformer
        for a, b in zip(list_smiles, list_compound_name):
            
            counter +=1
            output_filename = f"{output_folder}/biotransformer_output_{counter}.csv"

            try:
                os.remove("biotransformer_temp_output.csv")
            except FileNotFoundError:
                pass  # File does not exist, no need to remove
            
            try:
                if mode == 'standard' or 'btType':
                    biotransformcall = 'java -jar BioTransformer3.0_20230525.jar -k pred -b ' + type_of_biotransformation +' -ismi "' + a +'" -ocsv '+str(output_filename)+' -s '+str(number_of_steps)
                elif mode.str.startswith('-k'):
                    biotransformcall = 'java -jar BioTransformer3.0_20230525.jar -ismi "' + a +'" -ocsv '+str(output_filename)+' '+mode

                #print(biotransformcall)
                print('RUNNING VIRTUAL METABOLISATION WITH BIOTRANSFORMER')
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
                print('                    '+b)
                print('                    For:  '+a+' . Check the SMILES on http://biotransformer.ca and/or the GNPS library entry.')
                pass

        #Create a consensus name
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

        file_name = f"{output_name}_results_vm-NAP_BioTransformer"
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
        run_biotransformer3.markdown_link_biotransf = f"[View/Download the vm-NAP results from {output_name}](./{run_biotransformer3.file_name_biotransf})."
        run_biotransformer3.markdown_link_biotransf_nap = f"[View/Download the vm-NAP BioTransformer results for Network Annotation Propagation from {output_name}](./{run_biotransformer3.file_name_biotransf_nap})."
        run_biotransformer3.markdown_link_biotransf_sirius = f"[View/Download the vm-NAP BioTransformer results for SIRIUS from {output_name}](./{run_biotransformer3.file_name_biotransf_sirius})."

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e.output.decode()}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")