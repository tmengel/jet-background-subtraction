#!/bin/bash
rm slurm-batch/slurm-out/*
rm slurm-batch/slurm-err/*

export DATADIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets
export OUTDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/background-subtraction/classical-subtract/multiplicity-configs

for energy in 200GeV 2760GeV; do
    for radius in R02 R04 R06; do
        export INFILE=${DATADIR}/${energy}/${radius}/${energy}_Match_${radius}.root
        export OUTFILE=${OUTDIR}/${energy}_${radius}.root
        sbatch slurm-batch/config-job.sh $INFILE $OUTFILE;
    done
done