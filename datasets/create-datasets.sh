#!/bin/bash
make clean && make
#Clear slurm output files
rm $PWD/slurm-batch/slurm-err/*
rm $PWD/slurm-batch/slurm-out/*

export IN_DIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/jet-trees/root-files
export OUT_DIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/JetMult

# loop over energy directories (200GeV, 2760GeV) and jet radii (R02, R03, R04)
# for energy in 200GeV 2760GeV; do
#  for radius in R02 R04 R06; do 
#  sbatch slurm-batch/dataset-job.sh $IN_DIR/$energy/$radius $OUT_DIR
#  done
# done

sbatch slurm-batch/dataset-job.sh $IN_DIR/200GeV/R04 $OUT_DIR