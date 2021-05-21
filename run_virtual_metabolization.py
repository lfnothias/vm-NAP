from subprocess import call
import pandas as pd
import numpy as np
import sygma
from rdkit import Chem
import os
import sys

# Run SyGMa on two lists of SMILES and compound name
def run_sygma(list_smiles, list_compound_name, phase_1_cycle, phase_2_cycle, top_sygma_candidates):
    
    scenario = sygma.Scenario([
    [sygma.ruleset['phase1'], int(phase_1_cycle)],
    [sygma.ruleset['phase2'], int(phase_2_cycle)]])
    
    df = pd.DataFrame(data=None)
    df2 = pd.DataFrame(data=None)

    for a, b in zip(list_smiles, list_compound_name):
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
            #print(df)
            #print(pathway)
            df['pathway'] = pathway
            df['Compound_Name'] = b
            df2 = df2.append(df[:top_sygma_candidates], ignore_index=True)

        except:
            raise

    # Create new column name
    df2['score'] = df2['score'].astype(str)
    df2['Compound_Name_SyGMa'] = df2['pathway'] + df2['score']+ '; '+ df2['Compound_Name'] 
    df2["Compound_Name_SyGMa"] = df2["Compound_Name_SyGMa"].str.replace("\n", "")

    run_sygma.df2 = df2
    
    print('Number of SyGMA candidates = '+str(df2.shape[0]))
    print('Number of unique SyGMA candidates = '+str(len(df2.metabolite.unique())))

# Run BioTransformer on two lists of SMILES and compound name
def run_biotransformer(list_smiles, list_compound_name, type_of_biotransformation, number_of_steps):
    
    # Prepare pandas tables
    df_bio = pd.DataFrame(data=None)
    df2_bio = pd.DataFrame(data=None)
    counter = 0
    print('======== Running BioTransformer takes approximatively 3-5 compounds per minute ========')
    print('     Number of compounds being computed =  '+str(len(list_smiles)))

    # Iterative into the lists and run BioTransformer
    for a, b in zip(list_smiles[:5], list_compound_name[:5]):
        
        biotransformcall = 'java -jar biotransformer-1.1.5.jar  -k pred -b '+ type_of_biotransformation +' -ismi ' + a +' -ocsv biotransformer_temp_output.csv -s '+str(number_of_steps)
        biotransformcall = biotransformcall.split() # because call takes a list of strings 

        try:
            call(biotransformcall)

            # Read the results
            df_bio = pd.read_csv('biotransformer_temp_output.csv', sep=',')
            # Add the parent compound name
            df_bio['Parent_Compound_Name'] = b 
            counter +=1
            print('Number of BioTransformer candidates for compound '+str(counter)+ ' = '+str(df_bio.shape[0]))
            # Append to the master table
            df2_bio = df2_bio.append(df_bio, ignore_index=True)

        except:
            raise

    #Create a consensus name
    df2_bio['Compound_Name_BioTransformer'] = df2_bio['Reaction'] + '; '+ df2_bio['Parent_Compound_Name'] + '; '+ df2_bio['Reaction ID'] + '; '+ df2_bio['Biosystem']
    run_biotransformer.df2_bio = df2_bio 
    
    print('Total number of BioTransformer candidates = '+str(df2_bio.shape[0]))
    
#Write out the results
def export_for_SIRIUS(df, string):
    df_new = df # new table
    # RENAME SMILES
    df_new.to_csv('results_virtual_metabo_'+string+'.tsv', sep = '\t', index = False)
    
    try:
        #SyGMa
        df_new = df_new[['metabolite', 'Compound_Name_SyGMa']]
        df_new = df_new.rename(columns={"metabolite": "Smiles", "Compound_Name_SyGMA": "name"})
    except:
        #BioTransformer
        df_new = df_new[['SMILES', 'Compound_Name_BioTransformer']]
        df_new = df_new.rename(columns={"SMILES": "Smiles", "Compound_Name_BioTransformer": "name"})
        
    df_new.to_csv('results_virtual_metabo_formatted_'+string+'_SIRIUS.tsv', sep = '\t', index = False)
    
def export_for_NAP(df,string):
    df_new = df # new table
    # RENAME SMILES
    df_new.to_csv('results_virtual_metabo_'+str(string)+'.tsv', sep = '\t', index = False)
    try:
        #SyGMa
        df_new = df_new[['metabolite', 'Compound_Name_SyGMa']]
    except:
        #BioTransformer
        df_new = df_new[['SMILES', 'Compound_Name_BioTransformer']]

    df_new.to_csv('results_virtual_metabo_formatted_'+string+'_NAP.tsv', sep = '\t', index = False, header= False)