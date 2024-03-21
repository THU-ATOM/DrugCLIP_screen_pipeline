# DrugCLIP_screen_pipeline
A hybrid pipeline to screen compounds with DrugCLIP and Schrodinger

### 0. Requirements
First, you need [the Schrodinger suite](https://newsite.schrodinger.com/) for docking and other preprocessing.

Then, prepare a python environment with:

```
pip install pandas numpy biopython lmdb rdkit
```

You also need [the DrugCLIP model](https://github.com/bowen-gao/DrugCLIP), please see the other repo. However, you can integrate any virtual screening model into our pipeline.

### 1. Prepare chemical library for DrugCLIP

SDF files are usually provided by chemical suppliers like [ChemDiv](https://www.chemdiv.com/catalog/diversity-libraries/), [Enamine](https://enamine.net/compound-libraries/diversity-libraries) or [LifeChemicals](https://lifechemicals.com/screening-libraries/pre-plated-diversity-sets).

To convert SDF files into LMDB files that DrugCLIP model can process, you can put all your files into one folder and run:

```
python SDF2lmdb.py your_sdf_folder your_lmdb_file
```

### 2. Prepare protein (holo) pockets for DrugCLIP

DrugCLIP needs pre-defined pockets for screening. We recommend you use experimentally solved ligands to define the pocket, even if you would like to try DrugCLIP with AlphaFold2 models (you can align models first). Tools like [Fpocket](https://github.com/Discngine/fpocket) can be used, but only half of those pockets are precise enough for DrugCLIP.

After downloading your receptor-ligand structures from the PDB database, you should rename them as **PDBName_HetID.pdb**, where HetID is the molecule that defines the pocket. Put all PDB files into one folder, and run:

```
python pocket_from_pdb.py your_pdb_folder your_lmdb_name
```

The output LMDB file is located in the same folder as your PDB files

### 3. Virtual screening with DrugCLIP

You can use any other tools for virtual screening, but for further processing, you need to store your results as:

```
MolName,ChemSupplier,Score,SMILES
```

### 4. The novelty filter and clustering

We usually want molecules with novel structures (or cores), and for wet-lab screening, we often cannot afford a lot of molecules with similar structures. To remove molecules that are similar to known binders, you need to download activity data from the [ChEMBL](https://www.ebi.ac.uk/chembl/) website. Together with clustering, run:

```
python filter_cluster_pick.py screen_results ChEMBL_folder output_dir
```

You can also change other parameters like the fingerprint type for the novelty filter and clustering, or similarity cutoffs. For details, you might need to read the code.

### 5. Docking and scoring

As an AI-empowered virtual screening method, DrugCLIP can be hacked by out-of-distribution (OOD) molecules. Therefore, we use molecular docking as a physical-driven verification of screening results. In this step, only several hundred molecules are docked, so it is usually finished within an hour.

```
python glide_docking.py your_pdb_folder clustered_mols docking_outputs summerized_outputs
```

If your targets contain non-protein critical components in the pocket, like calcium cations, you need to preserve them in the docking grid by setting the **--keepHet** argument.

As for parallel computation, we usually recommend only one process for ligand preparation, but you might need to change it if you want to deal with a larger list. If you have more than one grid file (PDB target) for docking, please make sure:

```
Physical_Core_Number >= min(max_process,pdb_num) * ncpu_pergrid
```

### 6. Recycle to increase pocket diversity

DrugCLIP is designed to handle and ensemble multiple pockets at the same time. If you have only one holo PDB structure and your ligand cannot fully occupy the cavity, you can use this recycling strategy. Extract new pockets with the docking results from the previous step, and do **step 3-5** again:

```
#extract pockets with docking results
python pocket_from_sdf.py cleaned_pdb sdf_file your_lmdb_file
```

### 7. What makes a good candidate?

We recommend molecules with DrugCLIP zscore larger than 3, and docking score smaller than -6.
