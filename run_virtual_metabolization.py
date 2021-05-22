from subprocess import call
import pandas as pd
import numpy as np
import sygma
from rdkit import Chem
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)
import os
import sys

# Run SyGMa on two lists of SMILES and compound name in batch
def run_sygma_batch(list_smiles, list_compound_name, phase_1_cycle, phase_2_cycle, top_sygma_candidates, output_name, batch_size):
    
    print('=== Starting SyGMa computation ===')
    print('Number of compounds = '+str(len(list_smiles)))
    print('Batch_size = '+str(batch_size))
    print('If you are running many compounds or cycles, and maxing out RAM memory available, you can decreased the batch size. Otherwise the value can be increased for faster computation.')
    print('======')
    print('Please wait')
    print('======')
          
    scenario = sygma.Scenario([
    [sygma.ruleset['phase1'], int(phase_1_cycle)],
    [sygma.ruleset['phase2'], int(phase_2_cycle)]])
    
    df = pd.DataFrame(data=None)
    df2 = pd.DataFrame(data=None)
    df.to_csv('sygma_temp_output.csv', sep='\t')

    batch_counter = 0

    for i in range(0, len(list_smiles), batch_size):
        list_smiles_batch =  list_smiles[i:i + batch_size]
        list_compound_name_batch = list_compound_name[i:i + batch_size]

        for a, b in zip(list_smiles_batch, list_compound_name_batch):
            try:
                pathway = []
                mol = Chem.MolFromSmiles(a)
                metabolic_tree = scenario.run(mol)

                #Get the score
                metabolic_tree.calc_scores()

                #Get the pathway
                metabolites = metabolic_tree.to_smiles()

                # Get the metabolic pathway
                metabolite_info = metabolic_tree.to_list()
                for k in metabolite_info[1:]:
                    #print(k['SyGMa_pathway'])
                    pathway.append(k['SyGMa_pathway'])

                #Summarize results
                df = pd.DataFrame(metabolites[1:],columns=metabolites[0])
                df['parent'] = (metabolites[0][0])
                df.columns.values[0] = 'metabolite'
                df.columns.values[1] = 'score'
                df['score'] = df['score'].round(3)
                df['pathway'] = pathway
                df['pathway'] = df['pathway'].str[:75]
                df['Compound_Name'] = b[:50]
                df2 = df2.append(df[:top_sygma_candidates], ignore_index=True)

            except:
                raise
                
        df2.to_csv('sygma_temp_output.csv', sep='\t', index=False)
        df = pd.read_csv('sygma_temp_output.csv', sep='\t')
        
        batch_counter += 1
        print('Batch '+str(batch_counter)+' was compeleted')
    
    # Create new column name
    df = pd.DataFrame(data=None)
    df2 = pd.read_csv('sygma_temp_output.csv', sep='\t')
    df2['score'] = df2['score'].astype(str)
    df2['Compound_Name_SyGMa'] = df2['pathway'] + df2['score']+ '; '+ df2['Compound_Name'] 
    df2["Compound_Name_SyGMa"] = df2["Compound_Name_SyGMa"].str.replace("\n", "")

    print('Number of SyGMA candidates = '+str(df2.shape[0]))
    print('Number of unique SyGMA candidates = '+str(len(df2.metabolite.unique())))
    print('===== COMPLETED =====')
    df2.to_csv(output_name, sep='\t', index = False)
    
    df2 = pd.DataFrame(data=None)


# Run BioTransformer on two lists of SMILES and compound name
def run_biotransformer(list_smiles, list_compound_name, type_of_biotransformation, number_of_steps, output_name):
    
    # Prepare pandas tables
    df_bio = pd.DataFrame(data=None)
    df2_bio = pd.DataFrame(data=None)
    counter = 0
    print('======== Running BioTransformer takes approximatively 3-5 compounds per minute ========')
    print('     Number of compounds being computed =  '+str(len(list_smiles)))

    # Iterative into the lists and run BioTransformer
    for a, b in zip(list_smiles, list_compound_name):
        
        biotransformcall = 'java -jar biotransformer-1.1.5.jar  -k pred -b '+ type_of_biotransformation +' -ismi ' + a +' -ocsv biotransformer_temp_output.csv -s '+str(number_of_steps)
        biotransformcall = biotransformcall.split() # because call takes a list of strings 

        try:
            call(biotransformcall)

            # Read the results
            df_bio = pd.read_csv('biotransformer_temp_output.csv', sep=',')
            # Add the parent compound name
            df_bio['Parent_Compound_Name'] = b 
            counter +=1
            print('          BioTransformer candidates for compound n'+str(counter)+ ' = '+str(df_bio.shape[0]))
            # Append to the master table
            df2_bio = df2_bio.append(df_bio, ignore_index=True)

        except:
            raise

    #Create a consensus name
    df2_bio['Compound_Name_BioTransformer'] = df2_bio['Reaction'] + '; '+ df2_bio['Parent_Compound_Name'] + '; '+ df2_bio['Reaction ID'] + '; '+ df2_bio['Biosystem']
    
    print('===== COMPLETED =====')
    print('Total number of BioTransformer candidates = '+str(df2_bio.shape[0]))
    df2_bio.to_csv(output_name, sep='\t', index = False)
    
#Write out the results
def export_for_SIRIUS(input_df_name):
    df_new = pd.read_csv(input_df_name, sep='\t')
    
    try:
        #SyGMa
        df_new = df_new[['metabolite', 'Compound_Name_SyGMa']]
        df_new = df_new.rename(columns={"metabolite": "Smiles", "Compound_Name_SyGMa": "name"})
    except:
        #BioTransformer
        df_new = df_new[['SMILES', 'Compound_Name_BioTransformer']]
        df_new = df_new.rename(columns={"SMILES": "Smiles", "Compound_Name_BioTransformer": "name"})
        
    df_new.to_csv(input_df_name[:-4]+'_SIRIUS.tsv', sep = '\t', index = False)
    
def export_for_NAP(input_df_name):
    df_new = pd.read_csv(input_df_name, sep='\t')

    try:
        #SyGMa
        df_new = df_new[['metabolite', 'Compound_Name_SyGMa']]
    except:
        #BioTransformer
        df_new = df_new[['SMILES', 'Compound_Name_BioTransformer']]

    df_new.to_csv(input_df_name[:-4]+'_NAP.tsv', sep = '\t', index = False, header= False)