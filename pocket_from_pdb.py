import os
from Bio.PDB import PDBParser,Chain,Model,Structure
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import is_aa
from Bio.PDB.Residue import DisorderedResidue,Residue
from Bio.PDB.Atom import DisorderedAtom
import warnings
from Bio.PDB.StructureBuilder import PDBConstructionWarning
from tqdm import tqdm
import numpy as np
import lmdb
import numpy as np
import pickle
import re

def write_lmdb(data, lmdb_path, num):
    #resume
    
    env = lmdb.open(lmdb_path, subdir=False, readonly=False, lock=False, readahead=False, meminit=False, map_size=1099511627776)
    with env.begin(write=True) as txn:
        for d in data:
            txn.put(str(num).encode('ascii'), pickle.dumps(d))
            num += 1
    
    return num

warnings.filterwarnings(
    action='ignore',
    category=PDBConstructionWarning)
def extract_lig_recpt(biopy_model,ligname):
    lig_list = []
    tmp_chain = Chain.Chain('A') 
    rid = 0
    for chain in biopy_model:
        for res in chain:
            res.detach_parent()
            if not(is_aa(res,standard=True)):
                if res.resname == ligname:
                    lig_list.append((chain.id+'_'+str(res.id[1]),res))
                continue
            if res.is_disordered():
                if isinstance(res,DisorderedResidue):
                    res = res.selected_child
                    res.id = (res.id[0],rid,res.id[2])
                    rid += 1 
                    tmp_chain.add(res.copy())      
                else:
                    new_res = Residue(res.id,res.resname,res.segid)
                    for atom in res:
                        if isinstance(atom,DisorderedAtom):  
                            atom.selected_child.disordered_flag = 0
                            new_res.add(atom.selected_child.copy())
                        else:
                            new_res.add(atom)
                    res = new_res
                    res.id = (res.id[0],rid,res.id[2])
                    rid += 1 
                    tmp_chain.add(res.copy())       
            else:
                res.id = (res.id[0],rid,res.id[2])
                rid += 1 
                tmp_chain.add(res.copy())
    # tmp_structure = Structure.Structure('tmp')
    # tmp_model = Model.Model(0)
    # tmp_structure.add(tmp_model) 
    # tmp_model.add(tmp_chain)
  
    return tmp_chain,lig_list

def get_binding_pockets(biopy_chain,liglist):
    pockets = []
    for n,lig in liglist:
            lig_coord = np.array([i.coord for i in lig.get_atoms() if i.element!='H'])
            tmp_chain = Chain.Chain('A') 
            for res in biopy_chain:
                res_coord = np.array([i.get_coord() for i in res.get_atoms()])
                dist = np.linalg.norm(res_coord[:,None,:]-lig_coord[None,:,:],axis=-1).min()
                if dist<=6:
                    tmp_chain.add(res.copy())
            pockets.append((n,tmp_chain))
    return pockets  

def pocket2lmdb(pocket_name,biopy_chain,pdb_name):
    recpt = list(biopy_chain.get_atoms())
    pocket_atom_type = [x.element for x in recpt if x.element!='H']
    pocket_coord = [x.coord for x in recpt if x.element!='H']
    print(pdb_name+'_'+pocket_name,len(pocket_atom_type))
    return {
        'pocket': pdb_name+'_'+pocket_name,
        'pocket_atoms': pocket_atom_type,
        'pocket_coordinates': pocket_coord
    }


def process_one_pdbdir(dirs,name='pocket'):
    all_pocket = []
    for d in os.listdir(dirs):
        if '.pdb' in d and '_clean' not in d:
            try:
                p = PDBParser()
                model = p.get_structure('0',os.path.join(dirs,d))[0]  
                tmp_chain,liglist = extract_lig_recpt(model,re.search(r'_.*\.',d)[0][1:-1])
                pocket = get_binding_pockets(tmp_chain,liglist)
                pocket = [pocket2lmdb(n,p,d.split('.')[0]) for n,p in pocket]
                all_pocket += pocket
            except:
                pass
    write_lmdb(all_pocket,os.path.join(dirs,f'{name}.lmdb'),0)
if __name__ == '__main__':
    #args
    import argparse
    parser = argparse.ArgumentParser(description='extract pocket from pdb')
    parser.add_argument('dir', type=str, help='pdb dir with file names of name_hetid')
    parser.add_argument('--name', type=str, help='name of the output lmdb',default='pocket')
    args = parser.parse_args()

    process_one_pdbdir(args.dir,args.name)