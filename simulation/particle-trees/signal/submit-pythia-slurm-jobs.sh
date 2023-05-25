#!/bin/bash
echo "Submitting jobs to generate Pythia pp data"
echo $1
echo $2
for i in {0..24}; do sbatch $PWD/slurm-batch/pythia-job.sh $1 $2 $i; done