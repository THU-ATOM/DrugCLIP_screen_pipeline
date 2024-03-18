import multiprocessing
from rdkit import Chem
from tqdm import tqdm
import pickle
import lmdb
import numpy as np
import os
from rdkit.Chem import AllChem
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

def write_lmdb(data, lmdb_path, num):
    #resume
    env = lmdb.open(lmdb_path, subdir=False, readonly=False, lock=False, readahead=False, meminit=False, map_size=1099511627776)
    with env.begin(write=True) as txn:
        for d in data:
            txn.put(str(num).encode('ascii'), pickle.dumps(d))
            num += 1
    return num

def gen_conf(mol):
    try:
        # Generate a single conformation for the molecule
        mol = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(mol, numConfs=1,numThreads=1,pruneRmsThresh=1,maxAttempts=1000,useRandomCoords=False)
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1)
        mol = Chem.RemoveHs(mol)
        if mol.GetNumConformers() == 0:
            return None
        else:
            return {'coordinates': [np.array(mol.GetConformer(i).GetPositions()) for i in range(mol.GetNumConformers())], 'atoms': [a.GetSymbol() for a in mol.GetAtoms()], 'smi':Chem.MolToSmiles(mol), 'IDs':getID(mol)}
    except Exception as e:
        # Handle any errors that occur during processing
        print(f"Error processing molecule: {e}")
        return None

def getID(mol):
    for id_key in ['IDNUMBER','Catalog ID','ID']:
        try:
            return mol.GetProp(id_key)
        except:
            pass
    raise ValueError('No ID found in molecule')

def process_sdf_file(sdf_file,n_cpu=32):
    subset = os.path.basename(sdf_file).split('.')[0]
    # Read the SDF file and convert all molecules to RDKit mol objects
    suppl = Chem.SDMolSupplier(sdf_file)
    molecules = [mol for mol in suppl if mol is not None]

    #test if the sdf file can be identified
    try:
        getID(molecules[0])
    except:
        raise ValueError('No ID found in molecule')
    
    # Use multiple processes to process the molecules in parallel,use tqdm to show progress
    with multiprocessing.Pool(n_cpu) as pool:
        processed_molecules = list(tqdm(pool.imap(gen_conf, molecules), total=len(molecules)))

    # Discard errored molecules
    mol_data = [mol for mol in processed_molecules if mol is not None]
    [mol.update({'subset':subset}) for mol in mol_data]
    
    return mol_data

if __name__ == "__main__":
    #mol_data = process_sdf_file('/drug/DrugCLIP_chemdata_v2024/SDFiles/LC_10k_Pre_Plated_Diversity_Set_PS6.sdf')
    home = '/drug/DrugCLIP_chemdata_v2024/SDFiles'
    output = '/drug/DrugCLIP_chemdata_v2024/DrugCLIP_mols_v2024.lmdb'
    num = 0
    #list home, sort all files from smallest to largest
    filelist = os.listdir(home)
    filelist.sort(key=lambda x: os.path.getsize(os.path.join(home, x)))
    for f in filelist:
        mol_data = process_sdf_file(os.path.join(home,f))
        num = write_lmdb(mol_data,output, num)
        print(f'Finished processing {f}, {len(mol_data)} molecules added, currently {num} molecules in the LMDB.')
    