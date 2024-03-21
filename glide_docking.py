import os
from tqdm import tqdm
import warnings
import numpy as np
from Bio.PDB import PDBParser,Chain,Model,Structure
from Bio.PDB.PDBIO import PDBIO
import numpy as np
from Bio.PDB import is_aa
import warnings
from Bio.PDB.StructureBuilder import PDBConstructionWarning
from filter_cluster_pick import read_dcfile
import multiprocessing as mp
from functools import partial

warnings.filterwarnings(
    action='ignore',
    category=PDBConstructionWarning)

SCHRODINGER = '/opt/schrodinger2021-2'
MAE2SDF = os.path.join(os.path.dirname(os.path.abspath(__file__)),'mae2sdf.py')
print('using mae2sdf.py from',MAE2SDF)

def clean_protein(home,pdbfiles,keepHet=[]):
    recpt_name,hetid = pdbfiles.split('.')[0].split('_')
    f = open(os.path.join(home,pdbfiles))
    p = PDBParser()
    receptor_model = p.get_structure('recpt',f)[0]
    recpt_chain = Chain.Chain('A') 
    lig_coord = []
    rid = 1
    for chain in receptor_model:
        for res in chain:
            if res.resname == hetid:
                lig_coord.append(np.array([a.coord for a in res.get_atoms()]))
            if is_aa(res,standard=True) or res.resname in keepHet:
                res.detach_parent()
                res.id = (res.id[0],rid,res.id[2]) 
                recpt_chain.add(res)
                rid += 1
    if len(lig_coord) == 0:
        raise ValueError('Ligand not found')       
    tmp_structure = Structure.Structure('receptor')
    tmp_model = Model.Model(0)
    tmp_structure.add(tmp_model) 
    tmp_model.add(recpt_chain)  
    io = PDBIO()
    io.set_structure(tmp_structure)
    io.save(os.path.join(home,recpt_name+'_clean.pdb'))
    return lig_coord,recpt_name+'_clean.pdb'

def prepare_docking_grid(home,prot,lig_coord,overwrite=False):
    
    prot_fixed = os.path.join(home,os.path.splitext(prot)[0]+'_fixed.maegz')
    os.chdir(home)
    if os.path.exists(prot_fixed) and not overwrite:
        print('Using existed fixed pdb files!')
    else:
        #prepare protein
        prot_dir = os.path.join(home,prot)
        os.system(f'{SCHRODINGER}/utilities/prepwizard -j prepwizard_{os.path.splitext(prot)[0]} -watdist 5 -propka_pH 7.4 -rmsd 0.30 -HOST localhost:1 -NJOBS 1 -TMPLAUNCHDIR -ATTACHED -WAIT {prot_dir} {os.path.splitext(prot)[0]+"_fixed.maegz"}')
    grid_name = []
    for i,lig in enumerate(lig_coord):
        x,y,z = lig.mean(axis=0)
        grid_out = os.path.join(home,f'{os.path.splitext(prot)[0].split("_")[0]}_glidgrid{i}.zip')
        grid_name.append(grid_out)
        if os.path.exists(grid_out) and not overwrite:
            print('Using existed grid files!')
        else:
            #generate grid
            grid_conf = os.path.join(home,f'{os.path.splitext(prot)[0].split("_")[0]}_conf{i}.gridin')
            with open(grid_conf,'w') as f:
                f.write(
                    f'GRIDFILE {grid_out}\nRECEP_FILE  {prot_fixed}\nGRID_CENTER {x},{y},{z}'
                )
            os.system(f'{SCHRODINGER}/glide -WAIT -JOBNAME glide_grid_{os.path.splitext(prot)[0]} {grid_conf}')
    
    return grid_name

def prepare_ligand(dc_file,output_dir,overwrite=False,n_cpu=1):
    smile_string,subset_dict = read_dcfile(dc_file)
    with open(os.path.join(output_dir,'ligands.smi'),'w') as f:
        f.writelines(smile_string)
    os.chdir(output_dir)
    lig_conf = os.path.join(output_dir,'ligprep_conf.in')
    conformation_out = os.path.join(output_dir,'prepared_ligands.maegz')
    if os.path.exists(conformation_out) and not overwrite:
        print('Using existed prepared ligands!')
    else:
        with open(lig_conf,'w') as f:
            f.write(
                f'''INPUT_FILE_NAME {os.path.join(output_dir,'ligands.smi')}
OUT_MAE {conformation_out}
FORCE FIELD 16
EPIK yes
DETERMINE CHIRALITIES yes
RESPECT CHIRALITIES yes
IGNORE CHIRALITIES no
NUM STEREOISOMERS {32}
'''
            )
        os.system(f'{SCHRODINGER}/ligprep -inp {lig_conf} -HOST localhost:{n_cpu} -NJOBS 1 -JOBNAME ligprep -WAIT -LOCAL')


    return subset_dict 

def docking_with_grids(home,grid_list,PRECISION='SP',n_cpu=4,overwrite=False):
    glide_results = []
    os.chdir(home)
    for grid in grid_list:
        glide_out = os.path.join(home,f"gliddock_{os.path.basename(grid).split('.')[0]}_{PRECISION}_pv.maegz")
        glide_results.append(glide_out)
        if os.path.exists(glide_out) and not overwrite:
            print('Using existed docking results!')
        else:
            docking_conf = os.path.join(home,f"{os.path.basename(grid).split('.')[0]}_conf.dockingin")
            with open(docking_conf,'w') as f:
                f.write(
                    f'''GRIDFILE {grid}
PRECISION {PRECISION}
LIGANDFILE {os.path.join(home,'prepared_ligands.maegz')}
CALC_INPUT_RMS True
'''
                )    
            os.system(f"{SCHRODINGER}/glide -JOBNAME gliddock_{os.path.basename(grid).split('.')[0]}_{PRECISION} -HOST localhost:{n_cpu} -NJOBS {n_cpu} -WAIT -NOLOCAL {docking_conf}")    
            
    return glide_results

def docking_one_pair(protein_file,protein_folder,output_dir,keepHet=[],overwrite=False,PRECISION='SP',n_cpu=4):
    lig_coord,cleaned_name = clean_protein(protein_folder,protein_file,keepHet=keepHet)
    grid_name = prepare_docking_grid(protein_folder,cleaned_name,lig_coord,overwrite=overwrite)
    glide_results = docking_with_grids(output_dir,grid_name,PRECISION=PRECISION,n_cpu=n_cpu)
    return glide_results



def cross_docking(protein_folder,dcfile,output_dir,result_dir,keepHet=[],ncpu_ligand=1,ncpu_pergrid=4,max_process=16,overwrite=False,PRECISION='SP'):
    receptor_list =[]
    for f in os.listdir(protein_folder):
        if '.pdb' in f and '_clean' not in f:
            receptor_list.append(f)
    subset_dict = prepare_ligand(dcfile,output_dir,overwrite=overwrite,n_cpu=ncpu_ligand)
    max_process = min(max_process,len(receptor_list))
    #multiprocessing docing one pair of receptor_list, use tqdm
    with mp.Pool(max_process) as p:
        func = partial(docking_one_pair,protein_folder=protein_folder,output_dir=output_dir,keepHet=keepHet,overwrite=overwrite,n_cpu=ncpu_pergrid,PRECISION=PRECISION)
        glide_results = list(tqdm(p.imap(func,receptor_list),total=len(receptor_list)))
    with open(os.path.join(output_dir,'glide_results.txt'),'w') as f:
        for g in glide_results:
            for i in g:
                f.write(f'{i}\n')
    #remove all files from result_dir
    for f in os.listdir(result_dir):
        os.remove(os.path.join(result_dir,f))                   
    
    os.system(f"{SCHRODINGER}/run {MAE2SDF} {os.path.join(output_dir,'glide_results.txt')} {result_dir} 0")  
    g = open(os.path.join(result_dir,'final_result.csv'),'w')
    g.write('MolID,Subset,DCScore,GlideScore,BestDockingGrid,smiles\n')
    with open(os.path.join(result_dir,'ligands.aff'),'r') as f:
        for lines in f.readlines():
            MolID,gscore,bestgrid = lines.strip().split('\t')
            smiles,dcscore,subset = subset_dict[MolID]
            for s in subset:
                g.write(f'{MolID},{s},{dcscore},{gscore},{bestgrid},{smiles}\n')
    g.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Cross docking with Glide')
    parser.add_argument('protein_folder',type=str,help='Folder containing protein pdb files')
    parser.add_argument('dcfile',type=str,help='clustered molecules in drugclip file format')
    parser.add_argument('output_dir',type=str,help='Output directory for all docking results')
    parser.add_argument('result_dir',type=str,help='Output directory for summmarized results and sdf files')
    parser.add_argument('--keepHet',type=str,nargs='*',default=[],help='List of hetid to keep in the protein')
    parser.add_argument('--ncpu_ligand',type=int,default=1,help='Number of cpus for ligand preparation')
    parser.add_argument('--ncpu_pergrid',type=int,default=4,help='Number of cpus for each docking jobs')
    parser.add_argument('--max_process',type=int,default=16,help='Maximum number of docking task at the same time')
    parser.add_argument('--overwrite',action='store_true',help='Overwrite existing files')
    parser.add_argument('--PRECISION',type=str,default='SP',help='Glide precision')
    args = parser.parse_args()
    cross_docking(**vars(args))
