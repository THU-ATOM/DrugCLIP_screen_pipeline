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
from rdkit import Chem
from pocket_from_pdb import extract_lig_recpt,pocket2lmdb,write_lmdb

warnings.filterwarnings(
    action='ignore',
    category=PDBConstructionWarning)

def read_sdf(sdfile):
    #read a sdf file with rdkit and return all coord array in a list
    suppl = Chem.SDMolSupplier(sdfile)
    mol_coord = []
    for mol in suppl:
        if mol is not None:
            mol_coord.append(mol.GetConformer(0).GetPositions())
    return mol_coord

def get_binding_pockets(biopy_chain,coordlist):
    pockets = []
    for n,lig_coord in enumerate(coordlist):
            tmp_chain = Chain.Chain('A') 
            for res in biopy_chain:
                res_coord = np.array([i.get_coord() for i in res.get_atoms()])
                dist = np.linalg.norm(res_coord[:,None,:]-lig_coord[None,:,:],axis=-1).min()
                if dist<=6:
                    tmp_chain.add(res.copy())
            pockets.append((str(n),tmp_chain))
    return pockets  

def process_one_pair(pdbfile,sdfile,outputfile):
    p = PDBParser()
    model = p.get_structure('0',pdbfile)[0]  
    tmp_chain,_ = extract_lig_recpt(model,'nothing_here')    
    coordlist = read_sdf(sdfile)
    pocket = get_binding_pockets(tmp_chain,coordlist)
    pocket = [pocket2lmdb(n,p,os.path.basename(pdbfile).split('.')[0]) for n,p in pocket]
    write_lmdb(pocket,outputfile,0)

if __name__ == '__main__':
    process_one_pair('/newhome/jiayj/DrugCLIP_screen_pipeline/targets/GammaSecretase/Pocket2/E2012andL458_clean.pdb',
                     '/newhome/jiayj/DrugCLIP_screen_pipeline/targets/GammaSecretase/Pocket2_round1/outputs/E2012andL458_glidgrid0.sdf',
                     '/newhome/jiayj/DrugCLIP_screen_pipeline/targets/GammaSecretase/Pocket2/Pocket2_round2.lmdb')
    # import argparse
    # parser = argparse.ArgumentParser(description='extract binding pockets from pdb and sdf')
    # parser.add_argument('pdb', type=str, help='pdb file')
    # parser.add_argument('sdf', type=str, help='sdf file')
    # parser.add_argument('output', type=str, help='output lmdb')
    # args = parser.parse_args()
    # process_one_pair(args.pdb,args.sdf,args.output)
