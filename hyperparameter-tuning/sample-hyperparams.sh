#!/bin/bash
rm slurm-batch/slurm-out/*
rm slurm-batch/slurm-err/*


for i in {1..100}; do
    sbatch slurm-batch/sample-job.sh
done

