#!/bin/bash
make clean && make
#Clear slurm output files
rm $PWD/slurm-batch/slurm-err/*
rm $PWD/slurm-batch/slurm-out/*

export IN_DIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/jet-trees/root-files
export OUT_DIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets

# loop over energy directories (200GeV, 2760GeV) and jet radii (R02, R03, R04)
for energy in 200GeV 2760GeV; do
 for radius in R02 R04 R06; do 
 sbatch slurm-batch/dataset-job.sh $IN_DIR/$energy/$radius $OUT_DIR
 done
done
# old af

# export FILEDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees/event/root-files
# export OUTDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/root-files
# export RHIC_DATA=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees/event/root-files/200GeV
# export LHC_DATA=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees/event/root-files/2760GeV

# = not ready yet
## = complete


## sbatch slurm-batch/dataset-job.sh $RHIC_DATA/R02/AuAu-Formatted/200GeV_R02_AuAu_cent0to60_formatted.root AuAu_R02_unmatched $OUTDIR
# sbatch slurm-batch/dataset-job.sh $RHIC_DATA/R03/AuAu-Formatted/200GeV_R03_AuAu_cent0to60_formatted.root AuAu_R03_unmatched $OUTDIR
## sbatch slurm-batch/dataset-job.sh $RHIC_DATA/R04/AuAu-Formatted/200GeV_R04_AuAu_cent0to60_formatted.root AuAu_R04_unmatched $OUTDIR
# sbatch slurm-batch/dataset-job.sh $RHIC_DATA/R05/AuAu-Formatted/200GeV_R05_AuAu_cent0to60_formatted.root AuAu_R05_unmatched $OUTDIR
## sbatch slurm-batch/dataset-job.sh $RHIC_DATA/R06/AuAu-Formatted/200GeV_R06_AuAu_cent0to60_formatted.root AuAu_R06_unmatched $OUTDIR

## sbatch slurm-batch/dataset-job.sh $RHIC_DATA/R02/Matched-AuAu-Formatted/200GeV_R02_Matched_AuAu_cent0to60_formatted.root AuAu_R02 $OUTDIR
## sbatch slurm-batch/dataset-job.sh $RHIC_DATA/R03/Matched-AuAu-Formatted/200GeV_R03_Matched_AuAu_cent0to60_formatted.root AuAu_R03 $OUTDIR
## sbatch slurm-batch/dataset-job.sh $RHIC_DATA/R04/Matched-AuAu-Formatted/200GeV_R04_Matched_AuAu_cent0to60_formatted.root AuAu_R04 $OUTDIR
## sbatch slurm-batch/dataset-job.sh $RHIC_DATA/R05/Matched-AuAu-Formatted/200GeV_R05_Matched_AuAu_cent0to60_formatted.root AuAu_R05 $OUTDIR
## sbatch slurm-batch/dataset-job.sh $RHIC_DATA/R06/Matched-AuAu-Formatted/200GeV_R06_Matched_AuAu_cent0to60_formatted.root AuAu_R06 $OUTDIR

## sbatch slurm-batch/dataset-job.sh $LHC_DATA/R02/PbPb-Formatted/2760GeV_R02_PbPb_cent0to50_formatted.root PbPb_R02_unmatched $OUTDIR
## sbatch slurm-batch/dataset-job.sh $LHC_DATA/R03/PbPb-Formatted/2760GeV_R03_PbPb_cent0to50_formatted.root PbPb_R03_unmatched $OUTDIR
## sbatch slurm-batch/dataset-job.sh $LHC_DATA/R04/PbPb-Formatted/2760GeV_R04_PbPb_cent0to50_formatted.root PbPb_R04_unmatched $OUTDIR
# sbatch slurm-batch/dataset-job.sh $LHC_DATA/R05/PbPb-Formatted/2760GeV_R05_PbPb_cent0to50_formatted.root PbPb_R05_unmatched $OUTDIR
## sbatch slurm-batch/dataset-job.sh $LHC_DATA/R06/PbPb-Formatted/2760GeV_R06_PbPb_cent0to50_formatted.root PbPb_R06_unmatched $OUTDIR

## sbatch slurm-batch/dataset-job.sh $LHC_DATA/R02/Matched-PbPb-Formatted/2760GeV_R02_Matched_PbPb_cent0to50_formatted.root PbPb_R02 $OUTDIR
## sbatch slurm-batch/dataset-job.sh $LHC_DATA/R03/Matched-PbPb-Formatted/2760GeV_R03_Matched_PbPb_cent0to50_formatted.root PbPb_R03 $OUTDIR
## sbatch slurm-batch/dataset-job.sh $LHC_DATA/R04/Matched-PbPb-Formatted/2760GeV_R04_Matched_PbPb_cent0to50_formatted.root PbPb_R04 $OUTDIR
## sbatch slurm-batch/dataset-job.sh $LHC_DATA/R05/Matched-PbPb-Formatted/2760GeV_R05_Matched_PbPb_cent0to50_formatted.root PbPb_R05 $OUTDIR
## sbatch slurm-batch/dataset-job.sh $LHC_DATA/R06/Matched-PbPb-Formatted/2760GeV_R06_Matched_PbPb_cent0to50_formatted.root PbPb_R06 $OUTDIR


