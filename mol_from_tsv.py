import pandas as pd
import argparse
import multiprocessing as mp
import os
import pickle
from random import shuffle
import lmdb
import numpy as np
import pandas as pd
import rdkit
from tqdm import tqdm
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

ACTIVE_CODE = ['active','weak activity','slightly active','slight inhibition','potent inhibitor','partially active','partial antagonist','partial agonist','non-competitive antagonist','inverse agonist','irreversible antagonist','dose-dependent effect','antagonist','agonist','activator']
INACTIVE_CODE = ['no significant effect','no significant activity','no effect','no activity','no action','inhibition not detected','inactive']
def tsv2smiles(tsvfile,thres_active=10000,thres_inactive=20000):

    dataframe = pd.read_csv(tsvfile, sep="\t",quotechar = '"')
    ids = dataframe['Molecule ChEMBL ID'].tolist()
    smi = dataframe['Smiles'].tolist()
    relation = dataframe['Standard Relation'].tolist()
    value = dataframe['Standard Value'].tolist()
    unites = dataframe['Standard Units'].tolist()
    comments = dataframe['Comment'].tolist()
    active = set({})
    inactive = set({})
    not_determined= set({})
    for n,s,r,v,u,c in zip(ids,smi,relation,value,unites,comments):
        if str(c).lower() in ACTIVE_CODE:#respect expert opinions
            active.add((n,s))
        elif str(c).lower() in INACTIVE_CODE or 'not active' in str(c).lower() or 'no inhibit' in str(c).lower():
            inactive.add((n,s))
        else:        
            if u != 'nM':
                not_determined.add((n,s))
            else:
                if r != "'>'" and v <= thres_active:
                    active.add((n,s))
                elif r != "'<'" and v >= thres_inactive:
                    inactive.add((n,s))
                else:
                    not_determined.add((n,s))
    return active,inactive,not_determined

def write_lmdb(data, lmdb_path, num):
    #resume
    
    env = lmdb.open(lmdb_path, subdir=False, readonly=False, lock=False, readahead=False, meminit=False, map_size=1099511627776)
    with env.begin(write=True) as txn:
        for d in data:
            txn.put(str(num).encode('ascii'), pickle.dumps(d))
            num += 1
def gen_conformation(smiles, num_conf=1, num_worker=1):
    if len(smiles)>100:
        print("exceede max smiles lens", smiles)
        return None
    mol = Chem.MolFromSmiles(smiles)
    try:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(mol, numConfs=num_conf, numThreads=num_worker, pruneRmsThresh=1, maxAttempts=1000, useRandomCoords=False)
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=num_worker)
        mol = Chem.RemoveHs(mol)
    except:
        print("cannot gen conf", smiles)
        return None
    if mol.GetNumConformers() == 0:
        print("cannot gen conf", smiles)
        return None
    return mol

def convert_2Dmol_to_data(name,simles, num_conf=1, num_worker=1):
    #to 3D
    mol = gen_conformation(simles, num_conf, num_worker)
    if mol is None:
        return None
    coords = [np.array(mol.GetConformer(i).GetPositions()) for i in range(mol.GetNumConformers())]
    atom_types = [a.GetSymbol() for a in mol.GetAtoms()]
    return {'coordinates': coords, 'atoms': atom_types, 'smi': name+'_'+simles, 'mol': mol}

def multiprocess_confgen(smi_list,n_cpu=32):
    result = []
    tbar = tqdm(total=len(smi_list))
    def call_back(r):
        tbar.update(1)
        if r is not None:
            result.append(r)
    pool = mp.Pool(n_cpu)
    shuffle(smi_list)
    for n,s in smi_list:
        pool.apply_async(func=convert_2Dmol_to_data,args=(n,s,),callback=call_back)
    pool.close()
    pool.join()   
    return result 

def process_one_chembldir(dirs):
    active = set({})
    inactive = set({})
    not_determined= set({})
    for d in os.listdir(dirs):
        if '.tsv' in d:
            a,i,u = tsv2smiles(os.path.join(dirs,d))
            active.update(a)
            inactive.update(i)
            not_determined.update(u)
    active_x_inactive = active.intersection(inactive)
    active.difference_update(active_x_inactive)
    inactive.difference_update(active_x_inactive)
    not_determined.update(active_x_inactive)    
    active_mol = multiprocess_confgen(list(active))
    inactive_mol=multiprocess_confgen(list(inactive))
    not_determined_mol=multiprocess_confgen(list(not_determined))
    write_lmdb(active_mol, os.path.join(dirs,'active.lmdb'), 0)
    write_lmdb(inactive_mol, os.path.join(dirs,'inactive.lmdb'), 0)
    write_lmdb(not_determined_mol, os.path.join(dirs,'not_determined.lmdb'), 0)
    return len(active_mol),len(inactive_mol),len(not_determined_mol),len(active_x_inactive)

def PubChem_csv(csv_file):
    mols = set({})
    for lines in open(csv_file):
        cid = lines.split(',')[2]
        smi = lines.split(',')[3]
        if len(smi) != 0:
            mols.add((cid,smi))
    mol_data = multiprocess_confgen(list(mols))
    write_lmdb(mol_data, csv_file.replace('.csv','.lmdb'), 0)

if __name__ == '__main__':
    #arg 
    parser = argparse.ArgumentParser(description='convert tsv to 3D mol')
    parser.add_argument('dir', type=str, help='molecule information dir in ChEMBL format')
    args = parser.parse_args()
    process_one_chembldir(args.dir)


