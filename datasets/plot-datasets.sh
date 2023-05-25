#!/bin/bash
make clean && make
#Clear slurm output files
rm $PWD/slurm-batch/slurm-err/*
rm $PWD/slurm-batch/slurm-out/*

# export INDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets

# # loop over energy directories (200GeV, 2760GeV) and jet radii (R02, R03, R04)
# for energy in 200GeV 2760GeV; do
#     for radius in R02 R04 R06; do 
#         for filetype in Fake FullJets Match Missed PP; do
#             export DATADIR=${INDIR}/${energy}/${radius}/${energy}_${filetype}_${radius}.root
#             sbatch slurm-batch/plotting-job.sh $DATADIR;
#         done
#     done
# done

export INDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/results

# loop over energy directories (200GeV, 2760GeV) and jet radii (R02, R03, R04)
for energy in 200GeV 2760GeV; do
    for radius in R02 R04 R06; do 
        for filetype in Fake FullJets Match; do
            export DATADIR=${INDIR}/${energy}/${radius}/${energy}_${filetype}_${radius}_results.root
            sbatch slurm-batch/plotting-job.sh $DATADIR;
        done
    done
done
