import pandas as pd
import numpy as np
import math
import os
import sys


def print_compound_names(list_compounds):
    for item in sorted(list_compounds):
        print('\''+item+'\',')
    del list_compounds

    
    
def get_info_gnps_annotations(df_annotations, inchi_column, smiles_column, smiles_planar_column=False):
    #Get info on consolidated GNPS annotation table
    
    df_annotations_missing_structure = df_annotations[~df_annotations[inchi_column].str.startswith('InChI', na=False)]

    number_of_annotations = df_annotations.shape[0]
    number_of_annotations_without_structure = df_annotations_missing_structure.shape[0]

    print(str(number_of_annotations)+' annotations detected')
    list_unique_stereostructure = set(df_annotations[smiles_column])
    print('that corresponds to '+str(len(list_unique_stereostructure))+' unique stereostructures')
    
    if smiles_planar_column != False:
        list_unique_planarstructure = set(df_annotations[smiles_planar_column])
        print('that corresponds to '+(str(len(list_unique_planarstructure)))+' unique planar structures')

    print(' ==== WARNING =====')
    print(str(number_of_annotations_without_structure)+' annotations dont have a structure identifier and will be discarded from downstream processing, unless you do the following:')
    print('You can either update the GNPS library and rerun the GNPS job. Or you can provide a structure identifier in the dedicated cell below')    
    print('These are the compounds without structure identifiers:')
    list_missing_compounds = set(df_annotations_missing_structure['Compound_Name'])
    print_compound_names(list_missing_compounds)
    
    df_annotations = df_annotations[df_annotations[inchi_column].str.startswith('InChI', na=False)]

    return df_annotations


def print_compound_name_for_tags(df_annotations):
        mask = (df_annotations['tags'].astype('str').str.len() > 4)
        df = df_annotations.loc[mask]
        df = df.drop_duplicates(subset=['Compound_Name'])
        return df.sort_values(['tags'])[['tags','Compound_Name']]
    

def df_annotations_filtering(df_annotations, compound_name=False, tags=False):
    #If compound names or tags are available, we generate subtables
    if compound_name != False:
        df_annotations_name = df_annotations[df_annotations.Compound_Name.isin(compound_name)]
    if tags != False:
        df_annotations_tags = df_annotations[df_annotations.tags.isin(tags)]
     
    if compound_name != False and tags != False :
        df_annotations = pd.concat([df_annotations_name, df_annotations_tags], ignore_index=True)
        
    elif compound_name != False:
        df_annotations = df_annotations_name
        
    elif tags != False:
        df_annotations = df_annotations_tags
    else:
        print('No Compound_Name or Tags filter used')
        
    return df_annotations



def prepare_for_virtual_metabolization(df_annotations, compound_name, smiles_planar_column,  smiles_column=False, drop_duplicated_structure = True, use_planar_structure = True):
    
    # Input: consolidate GNPS annotation table
    # Do some duplicate filtering and select planar or stereo structure
    
    print('Number of spectral library annotations = '+str(df_annotations.shape[0]))
    
    if use_planar_structure == True:
        df_annotations = df_annotations[df_annotations[smiles_planar_column].str.contains('nan') == False]
        print('Number of spectral annotations with planar SMILES/InChI = '+str(df_annotations.shape[0]))

        if drop_duplicated_structure == True:
            try: 
                df_annotations = df_annotations.sort_values(by=['MQScore'], ascending=False)
            except:
                pass
                
            df_annotations = df_annotations.drop_duplicates(subset=smiles_planar_column, keep='first')

        list_compound_name = list(df_annotations[compound_name])
        list_smiles = list(df_annotations[smiles_planar_column])
        print('Number of unique planar SMILES considered = '+str(len(list_smiles)))
    
    else:
        df_annotations = df_annotations[df_annotations[smiles_column].str.contains('nan') == False]
        print('Number of spectral annotations with valid SMILES or InChI = '+str(df_annotations.shape[0]))
        
        if drop_duplicated_structure == True:
            try: 
                df_annotations = df_annotations.sort_values(by=['MQScore'], ascending=False)
            except:
                pass
            
            df_annotations = df_annotations.drop_duplicates(subset=smiles_column, keep='first')

        list_compound_name = list(df_annotations[compound_name])
        list_smiles = list(df_annotations[smiles_column])
        print('Number of unique SMILES = '+str(len(list_smiles)))
        
    prepare_for_virtual_metabolization.list_compound_name = list_compound_name
    prepare_for_virtual_metabolization.list_smiles = list_smiles
    
    return df_annotations



def append_to_list_if_not_present(base_list, extra_list):
    #Append items if not already in the list
    print('Initial number of items in the list = '+str(len(base_list)))
    for i in extra_list:
        if i not in base_list:
           base_list.append(i)
    print('Final number of items in the list = '+str(len(base_list)))
    
    return base_list