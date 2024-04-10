#!/bin/bash
rm slurm-batch/slurm-out/*
rm slurm-batch/slurm-err/*
# for i in 2 3 4; do
#     for j in A B C cos linear log; do
#         export INPUTFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v${i}/${j}_second_order_correlations_histos.root
#         for k in A B C cos linear log; do
#             export RESPONSEFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v${i}/${k}_second_order_correlations_histos.root
#             export OUTPUTFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v${i}/${j}_unfolded_with_${k}.root
#             # ROOT_COMMAND="root -l -b -q 'UnfoldFiles.C(\"$INPUTFILE\",\"$RESPONSEFILE\",\"$OUTPUTFILE\")'"
#             # echo $ROOT_COMMAND
#             sbatch slurm-batch/unfold-job.sh $INPUTFILE $RESPONSEFILE $OUTPUTFILE;
#         done
#     done
# done


# export INPUTFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v2/A_second_order_correlations_singlejet_areacut.root
# for j in A B C cos linear log; do
#     export RESPONSEFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v2/${j}_second_order_correlations_singlejet_areacut.root
#     export OUTPUTFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v2/A_unfolded_with_${j}_jetlevel_areacut_spectra.root
#     # sbatch slurm-batch/unfold-job.sh $INPUTFILE $RESPONSEFILE $OUTPUTFILE;
#     sbatch slurm-batch/job.sh $INPUTFILE $RESPONSEFILE $OUTPUTFILE;
# done


# export INPUTFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v2/A_second_order_correlations_singlejet_nocut.root
# for j in A B C cos linear log; do
#     export RESPONSEFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v2/${j}_second_order_correlations_singlejet_nocut.root
#     export OUTPUTFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v2/A_unfolded_with_${j}_jetlevel_nocut_spectra.root
#     # sbatch slurm-batch/unfold-job.sh $INPUTFILE $RESPONSEFILE $OUTPUTFILE;
#     sbatch slurm-batch/job.sh $INPUTFILE $RESPONSEFILE $OUTPUTFILE;
# done


# # for i in 2 3 4; do
# #     for j in A B C cos linear log; do
# #         export INPUTFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v${i}/${j}_second_order_correlations_singlejet_nocut.root
# #         for k in A B C cos linear log; do
# #             export RESPONSEFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v${i}/${k}_second_order_correlations_singlejet_nocut.root
# #             export OUTPUTFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v${i}/${j}_unfolded_with_${k}_jetlevel_nocut.root
# #             # ROOT_COMMAND="root -l -b -q 'UnfoldFiles.C(\"$INPUTFILE\",\"$RESPONSEFILE\",\"$OUTPUTFILE\")'"
# #             # echo $ROOT_COMMAND
# #             sbatch slurm-batch/unfold-job.sh $INPUTFILE $RESPONSEFILE $OUTPUTFILE;
# #         done
# #     done
# # done

# # for i in 2 3 4; do
# #     for j in A B C cos linear log; do
# #         export INPUTFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v${i}/${j}_second_order_correlations_singlejet_areacut.root
# #         for k in A B C cos linear log; do
# #             export RESPONSEFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v${i}/${k}_second_order_correlations_singlejet_areacut.root
# #             export OUTPUTFILE=/lustre/isaac/scratch/tmengel/JetVn/unfolding/root-files/twoparticle/v${i}/${j}_unfolded_with_${k}_jetlevel_areacut.root
# #             # ROOT_COMMAND="root -l -b -q 'UnfoldFiles.C(\"$INPUTFILE\",\"$RESPONSEFILE\",\"$OUTPUTFILE\")'"
# #             # echo $ROOT_COMMAND
# #             sbatch slurm-batch/unfold-job.sh $INPUTFILE $RESPONSEFILE $OUTPUTFILE;
# #         done
# #     done
# # done

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
    export OUTDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/unfold/results/AuAu_R0${i}/unfolding_output_new.root 
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
    export OUTDIR=/lustre/isaac/scratch/tmengel/jet-background-subtraction/unfold/results/PbPb_R0${i}/unfolding_output_new.root 
    # ROOT_COMMAND="root -l -b -q 'Unfold.C(\"$MATCHED\",\"$UNMATCHED\",\"$MISSED\",\"$FAKE\",\"$TRUTH\",\"$OUTDIR\")'"
    # eval $ROOT_COMMAND
    sbatch slurm-batch/unfold-job.sh $MATCHED $UNMATCHED $MISSED $FAKE $TRUTH $OUTDIR
done