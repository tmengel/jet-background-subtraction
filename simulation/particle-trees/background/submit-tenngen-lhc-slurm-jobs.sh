#!/bin/bash
echo "Submitting jobs to generate TennGen LHC data"
for i in {0..5}; do sbatch $PWD/slurm-batch/tenngen-lhc-job.sh $1 $i 0.9 1; done
