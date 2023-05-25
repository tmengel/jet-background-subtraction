#!/bin/bash
echo "Submitting jobs to generate TennGen data"
for i in {0,1,2,3}; do
#  for j in {0,1}; do 
sbatch $PWD/slurm-batch/tenngen-job.sh $1 $i 1.1 $2; done
# ; done