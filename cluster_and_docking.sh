#!/bin/bash

# Accepting the file path, directory, and name as command line arguments
file_path=$1 #drugclip output file
dir=$2  # directory for the whole target protein
name=$3 # name of current output files eg round2_out
chembl=$4 #directory name for chembl file eg ChEMBL
pdb=$5 #directory name for pdb file eg round2

# Creating the necessary directories
mkdir -p "$dir/$name/cluster"
mkdir -p "$dir/$name/docking"
mkdir -p "$dir/$name/output"

# Copying the file to dir/name/cluster/name.csv
cp "$file_path" "$dir/$name/cluster/$name.csv"

python filter_cluster_pick.py "$dir/$name/cluster/$name.csv" "$dir/$chembl" "$dir/$name/cluster"
wait

python glide_docking.py  "$dir/$pdb" "$dir/$name/cluster/dc_smiles_accepted.csv" "$dir/$name/docking" "$dir/$name/output"
wait