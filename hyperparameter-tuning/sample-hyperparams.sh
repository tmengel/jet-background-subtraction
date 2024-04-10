#!/bin/bash
rm slurm-batch/slurm-out/*
rm slurm-batch/slurm-err/*

cd slurm-batch
# for i in 1 10 20 30 40 50 60 70 80 90 100; do
#     for j in 1 10 20 30 40 50 60 70 80 90 100; do
#         for k in 1 10 20 30 40 50; do 
#         sbatch mse-job.sh $i $j $k
#         done
#     done
# done

for i in 1 10; do
    for j in 1 10; do
        for k in 1 10; do 
        sbatch mse-job.sh $i $j $k
        done
    done
done
