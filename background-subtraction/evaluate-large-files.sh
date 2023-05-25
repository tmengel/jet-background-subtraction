#!/bin/bash

rm slurm-batch/slurm-out/*
rm slurm-batch/slurm-err/*

export DATADIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets
export MODELDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/background-subtraction/ml-subtract/models
export RESULTDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/results

# export ENERGY=200GeV
# export SPECIES=AuAu
# for radius in R02 R04 R06; do
#     for filetype in Fake FullJets Match; do
#         export DNNMODEL=${MODELDIR}/${SPECIES}_${radius}_DNN_jet_pt_truth.root
#         export INFILE=${DATADIR}/${ENERGY}/${radius}/${ENERGY}_${filetype}_${radius}.root
#         export OUTFILE=${RESULTDIR}/${ENERGY}/${radius}/${ENERGY}_${filetype}_${radius}_results.root
#         sbatch slurm-batch/eval-job.sh $DNNMODEL $INFILE $OUTFILE;
#     done
# done

# export ENERGY=2760GeV
# export SPECIES=PbPb
# for radius in R02; do
#     for filetype in FullJets; do
#         export DNNMODEL=${MODELDIR}/${SPECIES}_${radius}_DNN_jet_pt_truth.root
#         export INFILE=${DATADIR}/${ENERGY}/${radius}/${ENERGY}_${filetype}_${radius}.root
#         export OUTFILE=${RESULTDIR}/${ENERGY}/${radius}/${ENERGY}_${filetype}_${radius}_results.root
#         sbatch slurm-batch/eval-job.sh $DNNMODEL $INFILE $OUTFILE;
#     done
# done

# export ENERGY=200GeV
# export SPECIES=AuAu
# for radius in R02 R04 R06; do
#     for filetype in Fake FullJets Match; do
#         export SNNMODEL=${MODELDIR}/${SPECIES}_${radius}_SNN_jet_pt_truth.root
#         export INFILE=${DATADIR}/${ENERGY}/${radius}/${ENERGY}_${filetype}_${radius}.root
#         export OUTFILE=${RESULTDIR}/${ENERGY}/${radius}/${ENERGY}_${filetype}_${radius}_results.root
#         sbatch slurm-batch/eval-job.sh $SNNMODEL $OUTFILE $OUTFILE;
#     done
# done

export ENERGY=2760GeV
export SPECIES=PbPb
for radius in R02; do
    for filetype in Fake FullJets Match; do
        export SNNMODEL=${MODELDIR}/${SPECIES}_${radius}_SNN_jet_pt_truth.root
        export INFILE=${DATADIR}/${ENERGY}/${radius}/${ENERGY}_${filetype}_${radius}.root
        export OUTFILE=${RESULTDIR}/${ENERGY}/${radius}/${ENERGY}_${filetype}_${radius}_results.root
        sbatch slurm-batch/eval-job.sh $SNNMODEL $OUTFILE $OUTFILE;
    done
done