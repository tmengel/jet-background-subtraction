#!/bin/bash
echo "Submitting jobs to generate Pythia pp data"
cd $PWD/signal
make clean
make
rm -rf slurm-batch/slurm-err/*
rm -rf slurm-batch/slurm-out/*

echo "Submitting jobs to generate Pythia pp data for RHIC"
./submit-pythia-slurm-jobs.sh 200 1000000
echo "Submitting jobs to generate Pythia pp data for LHC"
# ./submit-pythia-slurm-jobs.sh 2760 1000000

# cd ../background
# make clean
# make
# rm -rf slurm-batch/slurm-err/*
# rm -rf slurm-batch/slurm-out/*
# echo "Submitting jobs to generate TennGen data for RHIC"
# ./submit-tenngen-slurm-jobs.sh 1000000 1
# echo "Submitting jobs to generate TennGen data for LHC"
# ./submit-tenngen-lhc-slurm-jobs.sh 1000000

