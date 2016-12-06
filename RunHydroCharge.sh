#!/bin/bash
#purge all previous results
#What Ph to use?
pH=7.5

input=$1
results=$2
#TODO
#mkdir $results
#rm surfaces/*
#rm -rf surfacech/results/*


#you might want to comment the removal of contact maps below if you are re-running code -- calculation of these contact maps takes some time and it is better to cache them for re-running
#rm contact_maps/*

#[1/2] calculate surface charges for all pdbs found in abseqs
#cd surfacech
#python surfacer.py $pH $input
#mkdir -p ../$results/charge_per_residue
#cp results/*/*.tsv ../$results/charge_per_residue
#[2/2] calculate hydrophobicity stats for all pdbs -- store results in Results -- naming as described in the supplementary pdf where I explain how everything is calculated
#cd ../
python Annotator.py $pH $input $results

PML='
from glob import glob
for fname in glob(os.path.join(os.getcwd(), "annotated", "*.pdb")): cmd.load(fname)
cmd.zoom()
cmd.show_as("spheres")
cmd.spectrum("b")
'
for resdir in Results/*; do
  if [ -d "$resdir/annotated" ]; then
    echo "$PML" > $resdir/loadall.pml
  fi
done
