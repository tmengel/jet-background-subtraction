#!/bin/bash
rm slurm-batch/slurm-out/*
rm slurm-batch/slurm-err/*

for i in {2,4,6} ; do
    echo "Unfolding AuAu_R0${i}"
    export ENERGY=200GeV
    export DATADIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/${ENERGY}/R0${i}
    export RESULTDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/results/${ENERGY}/R0${i}
    export MATCHED=${RESULTDIR}/${ENERGY}_Match_R0${i}_results.root
    export UNMATCHED=${RESULTDIR}/${ENERGY}_FullJets_R0${i}_results.root
    export FAKE=${RESULTDIR}/${ENERGY}_Fake_R0${i}_results.root
    export MISSED=${DATADIR}/${ENERGY}_Missed_R0${i}.root
    export TRUTH=${DATADIR}/${ENERGY}_PP_R0${i}.root
    export OUTDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/unfold/results/AuAu_R0${i}/unfolding_output.root 
    # ROOT_COMMAND="root -l -b -q 'Unfold.C(\"$MATCHED\",\"$UNMATCHED\",\"$MISSED\",\"$FAKE\",\"$TRUTH\",\"$OUTDIR\")'"
    # eval $ROOT_COMMAND
    sbatch slurm-batch/unfold-job.sh $MATCHED $UNMATCHED $MISSED $FAKE $TRUTH $OUTDIR
done

for i in {2,4,6} ; do
    echo "Unfolding PbPb_R0${i}" 
    export ENERGY=2760GeV
    export DATADIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/${ENERGY}/R0${i}
    export RESULTDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/results/${ENERGY}/R0${i}
    export MATCHED=${RESULTDIR}/${ENERGY}_Match_R0${i}_results.root
    export UNMATCHED=${RESULTDIR}/${ENERGY}_FullJets_R0${i}_results.root
    export FAKE=${RESULTDIR}/${ENERGY}_Fake_R0${i}_results.root
    export MISSED=${DATADIR}/${ENERGY}_Missed_R0${i}.root
    export TRUTH=${DATADIR}/${ENERGY}_PP_R0${i}.root
    export OUTDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/unfold/results/PbPb_R0${i}/unfolding_output.root 
    # ROOT_COMMAND="root -l -b -q 'Unfold.C(\"$MATCHED\",\"$UNMATCHED\",\"$MISSED\",\"$FAKE\",\"$TRUTH\",\"$OUTDIR\")'"
    # eval $ROOT_COMMAND
    sbatch slurm-batch/unfold-job.sh $MATCHED $UNMATCHED $MISSED $FAKE $TRUTH $OUTDIR
done
