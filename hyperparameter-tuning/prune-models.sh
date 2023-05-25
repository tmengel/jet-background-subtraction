#!/bin/bash
rm slurm-batch/slurm-out/*
rm slurm-batch/slurm-err/*
rm plots/*

sbatch slurm-batch/prune-job.sh

