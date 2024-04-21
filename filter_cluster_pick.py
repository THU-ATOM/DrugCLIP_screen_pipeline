SCHRODINGER = '/opt/schrodinger2021-2'
import os
import pandas as pd

def extract_smi_from_dir(dirs,thres_active=100000):
    smi_list = set({})
    for f in os.listdir(dirs):
        if '.tsv' in f:
            dataframe = pd.read_csv(os.path.join(dirs,f), sep="\t",quotechar = '"')
            ids = dataframe['Molecule ChEMBL ID'].tolist()
            smi = dataframe['Smiles'].tolist()
            relation = dataframe['Standard Relation'].tolist()
            value = dataframe['Standard Value'].tolist()
            unites = dataframe['Standard Units'].tolist()
            comments = dataframe['Comment'].tolist()
            for n,s,r,v,u,c in zip(ids,smi,relation,value,unites,comments):
                if str(c).lower() in ['active','weak activity','slightly active','slight inhibition','potent inhibitor','partially active','partial antagonist','partial agonist','non-competitive antagonist','inverse agonist','irreversible antagonist','dose-dependent effect','antagonist','agonist','activator']: #respect expert opinions
                    smi_list.add(f'{s}\t{n}\n')
                else:        
                    if r != "'>'" and ((u == 'pM'and v <= thres_active*1000) or (u == 'nM'and v <= thres_active) or (u == 'uM'and v <= thres_active/1000)):
                        smi_list.add(f'{s}\t{n}\n')
    return list(smi_list)

def read_dcfile(dcfile):
    subset_dict = {}
    result = []
    with open(dcfile,'r') as f:
        for lines in f.readlines()[1:]: #with header
            ids,subset,zscore,smi=lines.strip().split(',')
            try:
                subset_dict[ids][-1].append(subset)
            except:
                subset_dict[ids] = (smi,f"{float(zscore):.2f}",[subset])
                result.append(f'{smi}\t{ids}\n')
    return result,subset_dict

def cluster_pick(home,smi,subset_dict,fp='maccs',dist=0.5,):
    os.chdir(home)
    os.system(f"{SCHRODINGER}/utilities/canvasFPGen -ismi {smi} -o {smi.replace('.smi','.fp')} -fptype {fp} -xp -compress")
    os.system(f"{SCHRODINGER}/utilities/canvasLC -ifp {smi.replace('.smi','.fp')} -ocsv LF_cluster.csv -group -dist {dist}")
    accepted = []
    with open(os.path.join(home,'LF_cluster.csv'),'r') as f:
        for lines in f.readlines():
            Name,Cluster,Leader = lines.strip().split(',')
            if Leader == '1':
                accepted.append(Name[1:-1])
    f = open(os.path.join(home,'dc_smiles_accepted.csv'),'w')
    f.write('ids,subset,zscore,smiles\n')
    for name in accepted:
        for subset in subset_dict[name][-1]:
            f.write(f'{name},{subset},{subset_dict[name][1]},{subset_dict[name][0]}\n')       
    f.close() 

def filter_cluster_pick(dcfile,exclude_folder,output_dir,fp_cluster='maccs',cluster_dict=0.5,fp_exclude='radial',cutoff=0.30,n_cpu=32):
    dc_mols,subset_dict = read_dcfile(dcfile)
    if os.path.isdir(exclude_folder):
        with open(os.path.join(output_dir,'dc_smiles.smi'),'w') as f:
            f.writelines(dc_mols)
        with open(os.path.join(output_dir,'active_smiles.smi'),'w') as f:
            f.writelines(extract_smi_from_dir(exclude_folder))
        os.chdir(output_dir)
        os.system(f"{SCHRODINGER}/utilities/canvasFPGen -JOB active_smiles -ismi {os.path.join(output_dir,'active_smiles.smi')} -o {os.path.join(output_dir,'active_smiles.fp')} -fptype {fp_exclude} -xp -compress -WAIT")
        os.system(f"{SCHRODINGER}/utilities/canvasFPGen -JOB dc_smiles -ismi {os.path.join(output_dir,'dc_smiles.smi')} -o {os.path.join(output_dir,'dc_smiles.fp')} -fptype {fp_exclude} -xp -compress -WAIT")    
        # os.system(f"{SCHRODINGER}/utilities/canvasDBCS -JOB  -HOST localhost:{n_cpu} -ifp {os.path.join(output_dir,'dc_smiles.fp')} -ifp2 {os.path.join(output_dir,'active_smiles.fp')} -n {n_mols} -d {cutoff} -o {os.path.join(output_dir,'DBCS.out')} -WAIT")   
        os.system(f"{SCHRODINGER}/utilities/canvasFPMatrix -JOB novelty_filter -HOST localhost:{n_cpu} -ifp {os.path.join(output_dir,'dc_smiles.fp')} -ifp2 {os.path.join(output_dir,'active_smiles.fp')} -filter {cutoff} -capcol 1 -ocsv {os.path.join(output_dir,'exclude_active.csv')} -WAIT")
        excluded = []
        with open(os.path.join(output_dir,'exclude_active.csv'),'r') as f:
            for lines in f:
                excluded.append(lines.split(',')[0])
    else:
        excluded = []
        Warning.warn(f'{exclude_folder} is not a directory, no exclusion will be performed.')
    f = open(os.path.join(output_dir,'dc_smiles_novel.smi'),'w') 
    g = open(os.path.join(output_dir,'dc_smiles_excluded.smi'),'w') 
    for lines in dc_mols:
        if lines.strip().split('\t')[1] in excluded:
            g.write(lines)
        else:
            f.write(lines)
    f.close()
    g.close()
    cluster_pick(
        home=output_dir,
        smi='dc_smiles_novel.smi',
        subset_dict = subset_dict,
        fp = fp_cluster,
        dist = cluster_dict,
    )
    print('finished') 


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='filter and cluster pick')
    parser.add_argument('dcfile', type=str, help='drugclip results')
    parser.add_argument('exclude_folder', type=str, help='folder with active molecules')
    parser.add_argument('output_dir', type=str, help='output directory')
    parser.add_argument('--fp_cluster', type=str, default='maccs', help='fingerprint for clustering')
    parser.add_argument('--cluster_dict', type=float, default=0.5, help='distance for clustering')
    parser.add_argument('--fp_exclude', type=str, default='radial', help='fingerprint for exclusion')
    parser.add_argument('--cutoff', type=float, default=0.30, help='cutoff for exclusion')
    parser.add_argument('--n_cpu', type=int, default=32, help='number of cpu')
    args = parser.parse_args()
    filter_cluster_pick(args.dcfile,args.exclude_folder,args.output_dir,args.fp_cluster,args.cluster_dict,args.fp_exclude,args.cutoff,args.n_cpu)

