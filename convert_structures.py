import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import SaltRemover
from molvs import standardize_smiles
from molvs import Standardizer
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)
import math

## Made by Louis Felix Nothias Feb 2021

def show_metadata_tools(table, metadata_prefix):
    #This prints the metadata columns per annotations using the respective prefix.
    metadata = []
    for x in table.columns:
        if x.startswith(metadata_prefix):
            metadata.append(str(x))
    print(metadata)

def consolidate_and_convert_structures(df, prefix, smiles, inchi=0):
	#This is the main funcition that takes a table, prefix for the column to be output, smiles name column, and inchi name column optionally
    if inchi == 0:
        df[smiles] = df[smiles].astype('str')

        print('Only SMILES were inputted')
        try: 
            make_mol_from_SMILES(df[smiles])
            mol_to_SMILES_iso(make_mol_from_SMILES.mol_list)
            mol_to_SMILES(make_mol_from_SMILES.mol_list)
            mol_to_INCHI(make_mol_from_SMILES.mol_list)
            mol_to_INCHIKEY(make_mol_from_SMILES.mol_list)

        except:
            raise
    else:
        print('Both SMILES and InChI were inputted')
        df[smiles] = df[smiles].astype('str')
        df[inchi] = df[inchi].astype('str')
        try:
            make_mol_from_SMILES(df[smiles])
            make_mol_from_INCHI(df[inchi])
            merge_mol(make_mol_from_SMILES.mol_list,make_mol_from_INCHI.mol_list)
            mol_to_SMILES_iso(merge_mol.consensus_mol_list)
            mol_to_SMILES(merge_mol.consensus_mol_list)
            mol_to_INCHI(merge_mol.consensus_mol_list)
            mol_to_INCHIKEY(merge_mol.consensus_mol_list)
            
        except: 
            raise

    df[str(prefix)+'Consol_SMILES_iso'] = mol_to_SMILES_iso.smiles_list
    df[str(prefix)+'Consol_SMILES'] = mol_to_SMILES.smiles_list
    df[str(prefix)+'Consol_InChIKey'] = mol_to_INCHIKEY.inchikey_list
    df[str(prefix)+'Consol_InChI'] = mol_to_INCHI.inchi_list
    
    mol_to_INCHI.inchi_list = []
    mol_to_INCHIKEY.inchikey_list = []
    mol_to_SMILES.smiles_list = []
    mol_to_SMILES_iso.smiles_list = []

    print('End')


def make_mol_from_SMILES(list_smiles):
	#Make a list of mol from SMILES list
    print('Converting SMILES to mol object')
    mol_list = []
    counter_success = 0
    counter_except = 0
    counter_not_available = 0

	# Run for loops for mol object from smiles list
    for x in list_smiles:
        if x is not np.nan and len(x) >5:
            try:
                mol_list.append(Chem.MolFromSmiles(remove_salt_from_SMILES(x)))
                #print('Conversion worked: '+str(x))
                counter_success += 1
            except:
                pass
                print('Conversion failed for: '+str(x))
                mol_list.append(np.nan)
                counter_except += 1
        else:
            #print('Not converted: '+str(x))
            mol_list.append(np.nan)
            counter_not_available += 1
            
    make_mol_from_SMILES.mol_list = mol_list
    mol_list = []
    print('Succesfully converted to mol object: '+str(counter_success))
    print('Exception to the parsing: '+str(counter_except))
    print('Not available: '+str(counter_not_available))


def make_mol_from_INCHI(list_inchi):
	#Make a list of mol from InChI list
    print('Converting INCHI to mol object')
    mol_list = []
    counter_success = 0
    counter_except = 0
    counter_not_available = 0
    
	# Run for loops for mol object from inchi list
    for x in list_inchi:
        if x is not np.nan and len(x) >7:
            try:
                mol_list.append(Chem.MolFromInchi(x))
                #print('Conversion worked: '+str(x))
                counter_success += 1
            except:
                print('Conversion failed for: '+str(x))
                mol_list.append(np.nan)
                counter_except += 1
                pass
        else:
            #print('Not converted: '+str(x))
            mol_list.append(np.nan)
            counter_not_available += 1
            
    make_mol_from_INCHI.mol_list = mol_list
    mol_list = []
    print('Succesfully converted to mol object: '+str(counter_success))
    print('Exception to the parsing: '+str(counter_except))
    print('Not available: '+str(counter_not_available))
    

def merge_mol(mol_list,mol_list2):
	#Consolidated a mol list from two lists of mols
    print('Consolidating the lists')
    consensus_mol_list = []
    counter_1 = 0
    counter_2 = 0
    counter_3 = 0
    for a, b in zip(mol_list, mol_list2):
        counter_3 +=1
        if a is not np.nan and a is not None:
            consensus_mol_list.append(neutralize_atoms_in_mol(mol_salt_remover(a)))
            counter_1 +=1
        elif b is not np.nan and b is not None:
            #b = mol_salt_remover(b)
            #b = neutralize_atoms_in_mol(b)
            consensus_mol_list.append(neutralize_atoms_in_mol(mol_salt_remover(b)))
            counter_2 +=1
        else:
            consensus_mol_list.append(np.nan)

    merge_mol.consensus_mol_list = consensus_mol_list
    consensus_mol_list = []
    print('Total mol object from the list 1 = '+str(counter_1))
    print('Mol object consolidated from list 2 = '+str(counter_2))
    print('Consolidated structures = '+str(counter_1+counter_2))


def mol_to_SMILES_iso(mol_list):
	#Take a list of mol objects and convert to canonical SMILES
	
    print('Converting mol objects to SMILES iso')
    smiles_list = []


    for x in mol_list:
        
        if x is not np.nan:
            try:
                smiles_list.append(standardize_smiles(remove_salt_from_SMILES(Chem.MolToSmiles(neutralize_atoms_in_mol(mol_salt_remover(x)), isomericSmiles=False))))
            except:

                smiles_list.append(np.nan)
                pass
        else:
            smiles_list.append(np.nan)
            pass

    mol_to_SMILES_iso.smiles_list = smiles_list
    smiles_list = []

def mol_to_SMILES(mol_list):
	#Take a list of mol objects and convert to SMILES with stereoconfiguration
    print('Converting mol objects to SMILES')
    smiles_list = []

    for x in mol_list:
        if x is not np.nan:
            try:
                smiles_list.append(standardize_smiles(remove_salt_from_SMILES(Chem.MolToSmiles(neutralize_atoms_in_mol(mol_salt_remover(x)), isomericSmiles=True))))
            except:
                pass
                smiles_list.append(np.nan)
        else:
            smiles_list.append(np.nan)
            pass
    
    mol_to_SMILES.smiles_list = smiles_list
    smiles_list = []

def mol_to_INCHIKEY(mol_list):
	#Take a list of mol objects and convert to INCHIKEY

    print('Converting mol objects to InChIKey')
    inchikey_list = []
    for x in mol_list:
        if x is not np.nan:
            try:
                x = mol_salt_remover(x)
                x = neutralize_atoms_in_mol(x)
                inchikey_list.append(Chem.MolToInchiKey(x))
            except:
                pass
                inchikey_list.append(np.nan)
        else:
            inchikey_list.append(np.nan)
            pass
        
    mol_to_INCHIKEY.inchikey_list = inchikey_list
    inchikey_list = []
    

def mol_to_INCHI(mol_list):
	#Take a list of mol objects and convert to INCHI
    print('Converting mol objects to InChI')
    inchi_list = []
    for x in mol_list:
        if x is not np.nan:
            try:
                x = mol_salt_remover(x)
                x = neutralize_atoms_in_mol(x)
                inchi_list.append(Chem.MolToInchi(x))
            except:
                pass
                inchi_list.append(np.nan)
        else:
            inchi_list.append(np.nan)
            pass
      
    mol_to_INCHI.inchi_list = inchi_list
    inchi_list = []  
    
#Remove salts from mol 
def mol_salt_remover(mol):
    remover = Chem.SaltRemover.SaltRemover(defnData="[Cl,F,Br,I,Na]")
    
    try:
        s = Standardizer()
        mol = s.standardize(mol)
        #return mol

    except:
        pass

    try:
        res, deleted = remover.StripMolWithDeleted(mol)
        neutralize_atoms(res)
        
        if len(deleted) > 0:
            print('  Salt(s) deleted in       : '+ str(Chem.MolToSmiles(mol)))
            print('  Remaining neutral residue: '+ str(Chem.MolToSmiles(res)))   
        res 
            
    except:
        pass

    return res

def neutralize_atoms_in_mol(mol):
    try: 
        pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()

    except:
        print('Neutralize atoms failed')
        print(mol)
        raise

    return mol

def remove_salt_from_SMILES(smiles):
    try:
        if '.' in smiles:
            print('  Salt(s) deleted in       : '+ str(smiles))
            smiles = smiles.split(".", 3)[0]
            smiles = standardize_smiles(smiles)
            mol = Chem.MolFromSmiles(smiles)
            neutralize_atoms_in_mol(mol)
            smiles = Chem.MolToSmiles(mol)
            print('  Remaining residue        : '+ str(smiles))
        return smiles

    except:
        raise
        
