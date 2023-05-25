#!/bin/bash

rm slurm-batch/slurm-out/*
rm slurm-batch/slurm-err/*

export DATADIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/background-subtraction/ml-subtract/results

for dataset in AuAu_R02 AuAu_R04 AuAu_R06 PbPb_R02 PbPb_R04 PbPb_R06; do
    for model in DNN SNN; do 
        export DATA=${DATADIR}/${dataset}_results.h5
        export OUTFILE=${dataset}_${model}
        sbatch slurm-batch/pysr-job.sh $DATA $OUTFILE $model;
    done
done
