import os
import multiprocessing as mp
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors
clean_up_list = []
ROS = '/drug/rosetta.binary.linux.release-315/main'
PY27 = '/home/jiayinjun/miniconda3/envs/py27/bin/python'
BCL = '/home/jiayinjun/flex_dock/bcl-4.3.1-Linux-x86_64/bcl.exe'
M2P = f'{ROS}/source/scripts/python/public/molfile_to_params.py'
CLP = f'{ROS}/tools/protein_tools/scripts/clean_pdb.py'
RSC = f'{ROS}/source/bin/rosetta_scripts.static.linuxgccrelease'
#get path of current file
ROP = os.path.join(os.path.dirname(os.path.abspath(__file__)),'RosettaLigandOptions.txt')
XML = os.path.join(os.path.dirname(os.path.abspath(__file__)),'flexible_docking.xml')

def rosetta_docking(lig_sdf_string,prot_name,work_dir,n_pose=5):
    os.chdir(work_dir)
    tmp_files = []
    lig_name = lig_sdf_string.split('\n')[0]
    try:
        with open(f'{lig_name}.sdf','w') as f:
            f.write(lig_sdf_string)
            f.write('$$$$\n')
        tmp_files.append(f'{lig_name}.sdf')
        #read sdf file with rdkit and calculate molecular weight
        mol = Chem.SDMolSupplier(f'{lig_name}.sdf')[0]
        mw = Descriptors.MolWt(mol)
        if mw>575 or mw<200: #docking relatively large mols to the pocket to enlarge it, avoid super large ligands
            raise ValueError('molecular weight out of range')
        os.system(f'''
{BCL} molecule:Filter -add_h -defined_atom_types \
  -3d -input_filenames {lig_name}.sdf \
  -output_matched {lig_name}.CLEANED.sdf \
  -output_unmatched {lig_name}.UNCLEANED.sdf -message_level Debug >/dev/null 2>&1
wait
{BCL} molecule:ConformerGenerator -rotamer_library cod \
  -top_models 100 -ensemble_filenames {lig_name}.CLEANED.sdf \
  -conformers_single_file {lig_name}.CLEANED.conf.sdf \
  -conformation_comparer 'Dihedral(method=Max)' 30 -max_iterations 1000 >/dev/null 2>&1
wait
{PY27} {M2P} -n LIG -p {lig_name} --chain=X --conformers-in-one-file {lig_name}.CLEANED.conf.sdf >/dev/null 2>&1
wait
cat {prot_name}_clean_A.pdb {lig_name}.pdb > {prot_name}_{lig_name}.pdb
'''
    )
        tmp_files.append(f'{lig_name}.CLEANED.sdf')
        tmp_files.append(f'{lig_name}.UNCLEANED.sdf')
        tmp_files.append(f'{lig_name}.CLEANED.conf.sdf')
        tmp_files.append(f'{lig_name}.params')
        tmp_files.append(f'{lig_name}.pdb')
        tmp_files.append(f'{lig_name}_conformers.pdb')
        tmp_files.append(f'{prot_name}_{lig_name}.pdb')
        #read the last line of {ligand_name}.params
        if os.path.exists(f'{lig_name}.params'):
            with open(f'{lig_name}.params','r') as f:
                last_line = f.readlines()[-1]
            assert last_line == f'PDB_ROTAMERS {lig_name}_conformers.pdb\n', 'ligand params not generated'
        else:
            raise ValueError('ligand params not generated')
        os.system(f"{RSC} @{ROP} -s {prot_name}_{lig_name}.pdb -extra_res_fa {lig_name}.params -nstruct {n_pose} -parser:protocol {XML} >/dev/null 2>&1")
        tmp_files.append(f'score.sc')
    except Exception as e:
        print(f'Error in docking {lig_name} to {prot_name}: {e}')
    for f in tmp_files:
        if os.path.exists(f):
            os.remove(f)    

def flexible_docking(protein_dir,ligand_dir,output_dir,n_pose=5,n_cpu=32):
    os.chdir(output_dir)
    prot_list=[]
    lig_list=[]
    clean_up_list=[]
    for f in os.listdir(ligand_dir):
        if f.endswith('.sdf'):
            protein_name = f.split('_')[0]
            os.system(f'{PY27} {CLP} {os.path.join(protein_dir,protein_name+"_clean.pdb")} A >/dev/null 2>&1')
            clean_up_list.append(os.path.join(output_dir,protein_name+'_clean_A.pdb'))
            clean_up_list.append(os.path.join(output_dir,protein_name+'_clean_A.fasta'))
            with open(os.path.join(ligand_dir,f),'r') as f:
                curr_lig = f.read().split('$$$$\n')[:-1]
            lig_list+=curr_lig
            prot_list+=[protein_name]*len(curr_lig)
    #run rosetta docking with multiprocessing, each process runs one ligand,
    #this is to avoid memory leak in rosetta
    pool = mp.Pool(n_cpu)
    tbar = tqdm(total=len(lig_list))
    def update(*a):
        tbar.update()
    for lig_sdf_string,prot_name in list(zip(lig_list,prot_list)):
        pool.apply_async(rosetta_docking,args=(lig_sdf_string,prot_name,output_dir,n_pose),callback=update)
    pool.close()
    pool.join()

    for f in clean_up_list:
        if os.path.exists(f):
            os.remove(f)
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='flexible docking with rosetta')
    parser.add_argument('protein_dir', type=str, help='protein dir')
    parser.add_argument('ligand_dir', type=str, help='ligand dir')
    parser.add_argument('output_dir', type=str, help='output dir')
    parser.add_argument('--n_pose', type=int, help='number of poses to generate for each ligand',default=5)
    parser.add_argument('--n_cpu', type=int, help='number of cpu to use',default=32)
    args = parser.parse_args()
    flexible_docking(args.protein_dir,args.ligand_dir,args.output_dir,args.n_pose,args.n_cpu)