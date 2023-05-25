#!/bin/bash
rm slurm-batch/slurm-out/*
rm slurm-batch/slurm-err/*

export DATADIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets
export OUTDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/background-subtraction/classical-subtract/multiplicity-configs

for energy in 200GeV 2760GeV; do
    for radius in R02 R04 R06; do
        for filetype in Fake FullJets Match; do
            export CONFIGFILE=${OUTDIR}/${energy}_${radius}.root
            export INFILE=${DATADIR}/${energy}/${radius}/${energy}_${filetype}_${radius}.root
            sbatch slurm-batch/subtract-job.sh $CONFIGFILE $INFILE;
        done
    done
done