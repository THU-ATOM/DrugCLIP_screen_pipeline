import os
from Bio.PDB.PDBIO import PDBIO
import numpy as np
from Bio.PDB import PDBParser,Chain,Atom,Residue

def pqr_parser(filename,return_score=False):
    with open(filename,'r') as f:
        data = f.readlines()
    coord = []
    for l in data:
        if "Pocket Score" in l:
            score = float(l.split()[-1])
        elif "Real volume" in l:
            volume = float(l.split()[-1])
        elif l[:4] == 'ATOM':
            coord.append(
                [float(l[30:38]),float(l[38:46]),float(l[46:54])]
            )
    coord = np.array(coord)
    if return_score:
        return coord,score,volume
    else:
        return coord


def build_complex(recpt_pdb,lig_coord):
    p = PDBParser()
    struct = p.get_structure('0',recpt_pdb) 
    new_chain = Chain.Chain('0')
    struct[0].add(new_chain)
    lig_res = Residue.Residue(('H_STP',1,' '), 'STP', ' ')
    new_chain.add(lig_res)
    for i,coord in enumerate(lig_coord):
        atom = Atom.Atom(f'F{i:03}',coord,1.0,1.0,' ',f'F{i:03}',i,'C')
        lig_res.add(atom)
    io = PDBIO()
    io.set_structure(struct)
    io.save(recpt_pdb.replace('.pdb','_STP.pdb'))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='build complex from pocket and receptor')
    parser.add_argument('receptor',help='receptor pdb file')
    parser.add_argument('ligand',help='ligand pqr file')
    args = parser.parse_args()
    lig_coord = pqr_parser(args.ligand)
    build_complex(args.receptor,lig_coord)    
