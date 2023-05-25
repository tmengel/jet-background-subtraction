#!/bin/bash
rm slurm-batch/slurm-out/*
rm slurm-batch/slurm-err/*

# create config files
python3 create-configs.py

for i in AuAu_R02 AuAu_R04 AuAu_R06 PbPb_R02 PbPb_R04 PbPb_R06; do
    sbatch slurm-batch/ml-job.sh Configs/$i.json;
done
