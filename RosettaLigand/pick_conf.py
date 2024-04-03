import numpy as np
import os
from tqdm import tqdm
from copy import deepcopy
import shutil
def read_pdb(pdb_file):
    #read all lines start from ATOM, extract x,y,z and add to a array, get score from the line start with Grid_score, and break
    with open(pdb_file,'r') as f:
        lines = f.readlines()
    xyz = []
    lig = []
    for line in lines:
        if line.startswith('ATOM') and not line.strip().endswith('H'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            xyz.append([x,y,z])
        elif line.startswith('HETATM') and not line.strip().endswith('H'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            lig.append([x,y,z])
        elif line.startswith('interface_delta_X'):
            score = float(line.split()[1])
            break
    return np.array(xyz),np.array(lig),score

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

def pick_conf(pdb_dir,output_dir,subname='',ref_pock=None,ppc_thres=4,interface_score_thres=-7,rmsd_thres=3.0,diff_num_thres=12):
    name = []
    coord = []
    score = []
    picked_pdb = []
    if ref_pock is not None:
        ref_coord = pqr_parser(ref_pock).mean(axis=0,keepdims=True)
    else:
        ref_coord = None
    for f in tqdm(os.listdir(pdb_dir)):
        if f.endswith('.pdb') and subname in f:
            xyz,lig_xyz,s = read_pdb(os.path.join(pdb_dir,f))
            if s>interface_score_thres:
                continue
            if ref_coord is not None and np.linalg.norm(lig_xyz-ref_coord,axis=-1).min()>ppc_thres:
                continue
            name.append(f)
            coord.append(xyz)
            score.append(s)
    coord = np.array(coord)
    score = np.array(score)
    name = np.array(name)
    while len(score)>0:
        #pick the first one
        idx = np.argmin(score)
        picked_pdb.append(deepcopy(name[idx]))
        distance = np.linalg.norm(coord-coord[idx],axis=-1)
        accepted_idx = ((distance>rmsd_thres).sum(axis=-1))>diff_num_thres
        coord = coord[accepted_idx]
        score = score[accepted_idx]
        name = name[accepted_idx]
    for f in picked_pdb:
        shutil.copy(os.path.join(pdb_dir,f),os.path.join(output_dir,f.replace('_','-').replace('.pdb','_LIG.pdb')))
    return picked_pdb
    
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='pick conformers from flexible docking')
    parser.add_argument('pdb_dir', type=str, help='pdb dir')
    parser.add_argument('output_dir', type=str, help='output dir')
    parser.add_argument('--pdb_prefix', type=str, help='only pdb files startwith the prefix will be used, default all pdb',default='')
    parser.add_argument('--ref_pock', type=str, help='reference fpocket .pqr file, used to exclude out-of-range docking pose',default=None)
    parser.add_argument('--ppc_thres', type=float, help='distance threshold to the reference pocket center, default 4A from pocketpicker paper',default=4)
    parser.add_argument('--interface_score_thres', type=float, help='interface score threshold, default -7',default=-7)
    parser.add_argument('--rmsd_thres', type=float, help='rmsd threshold, default 3.0',default=3.0)
    parser.add_argument('--diff_num_thres', type=int, help='minimum number of different atoms between chosen conformers, default 12',default=12)
    args = parser.parse_args()
    pick_conf(args.pdb_dir,args.output_dir,args.pdb_prefix,args.ref_pock,args.ppc_thres,args.interface_score_thres,args.rmsd_thres,args.diff_num_thres)
