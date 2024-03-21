import os
from schrodinger import structure,adapter
from rdkit import Chem

def process_multiple_mae(filelist,output_dir,gbsa):
    lig_dict = {}
    with open(filelist,'r') as f:
        maefiles = f.readlines()
    for f in maefiles:
        lig_dict = process_one_mae(f.strip(),lig_dict,gbsa)
    with open(os.path.join(output_dir,f'ligands.aff'),'w') as f:
        for name in lig_dict:
            f.write(f"{name}\t{lig_dict[name][0]}\t{lig_dict[name][2]}\n")
    
    for name in lig_dict:
        with open(os.path.join(output_dir,f'{lig_dict[name][2]}.sdf'),'a') as sdfFile:
            mol_block = Chem.MolToMolBlock(lig_dict[name][1])
            sdfFile.write(mol_block + '\n$$$$\n')

def process_one_mae(maefile,lig_dict,gbsa):
    if gbsa:
        e_name = 'r_psp_MMGBSA_dG_Bind'
    else:
        e_name = 'r_i_docking_score'
    grid_name = os.path.basename(maefile).split('_')[1:3]    
    grid_name = '_'.join(grid_name)
    indata =  structure.StructureReader(maefile)
    skip_first = 0
    for st in indata:   
        if skip_first:
            name = st.property['s_m_title']
            score = st.property[e_name]
            mol = adapter.to_rdkit(st)
            mol.SetProp('_Name',name)
            if name in lig_dict:
                if score<lig_dict[name][0]:
                    lig_dict.update({name:(score,mol,grid_name)})
                else:
                    pass
            else:
                lig_dict.update({name:(score,mol,grid_name)})
        else:
            skip_first = 1
    return lig_dict



if __name__ == '__main__':
    import sys
    process_multiple_mae(sys.argv[1],sys.argv[2],bool(int(sys.argv[3])))