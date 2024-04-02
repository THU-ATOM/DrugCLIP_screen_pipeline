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
    for line in lines:
        if line.startswith('ATOM') and not line.strip().endswith('H'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            xyz.append([x,y,z])
        elif line.startswith('Grid_score'):
            score = float(line.split()[1])
            break
    return np.array(xyz), score

def pick_conf(pdb_dir,output_dir,subname='',thres=(3,12)):
    name = []
    coord = []
    score = []
    picked_pdb = []
    for f in tqdm(os.listdir(pdb_dir)):
        if f.endswith('.pdb') and subname in f:
            xyz,s = read_pdb(os.path.join(pdb_dir,f))
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
        accepted_idx = ((distance>thres[0]).sum(axis=-1))>thres[1]
        coord = coord[accepted_idx]
        score = score[accepted_idx]
        name = name[accepted_idx]
    for f in picked_pdb:
        shutil.copy(os.path.join(pdb_dir,f),os.path.join(output_dir,f.replace('_0001.pdb','').replace('_','-')+'_LIG.pdb'))
    return picked_pdb
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='pick conformers from flexible docking')
    parser.add_argument('pdb_dir', type=str, help='pdb dir')
    parser.add_argument('output_dir', type=str, help='output dir')
    parser.add_argument('--pdb_prefix', type=str, help='only pdb files startwith the prefix will be used, default all pdb',default='')
    parser.add_argument('--threshold', type=float, nargs=2, help='threshold for picking conformers, (r,n) means at least n atoms have rmsd of r, default (3,12)',default=(3,12))
    pick_conf(parser.pdb_dir,parser.output_dir,parser.pdb_prefix,parser.threshold)
