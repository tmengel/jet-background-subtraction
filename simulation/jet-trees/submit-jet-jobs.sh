# rm slurm-batch/slurm-err/*
# rm slurm-batch/slurm-out/*
make clean
make
echo "Submitting jobs to jetfind and jetmatch."   
# modes

# for i in {0..24}; do
    # sbatch slurm-batch/jet-job.sh 200 $i 0.2 PP -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.2 FullJets -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.2 Match -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.2 Missed -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.2 Fake -1 1;
# done 

# for i in {0..24}; do
    # sbatch slurm-batch/jet-job.sh 200 $i 0.4 PP -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.4 FullJets -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.4 Match -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.4 Missed -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.4 Fake -1 1;
# done 

# for i in {0..24}; do
    # sbatch slurm-batch/jet-job.sh 200 $i 0.6 PP -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.6 FullJets -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.6 Match -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.6 Missed -1 1;
    # sbatch slurm-batch/jet-job.sh 200 $i 0.6 Fake -1 1;
# done 

# for i in {0..24}; do
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.2 PP -1 1;
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.2 FullJets -1 1;
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.2 Match -1 1;
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.2 Missed -1 1;
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.2 Fake -1 1;
# done 

# for i in {0..24}; do
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.4 PP -1 1;
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.4 FullJets -1 1;
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.4 Match -1 1;
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.4 Missed -1 1;
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.4 Fake -1 1;
# done 

for i in {0..24}; do
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.6 PP -1 1;
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.6 FullJets -1 1;
    # sbatch slurm-batch/jet-job.sh 2760 $i 0.6 Match -1 1;
    sbatch slurm-batch/jet-job.sh 2760 $i 0.6 Missed -1 1;
    sbatch slurm-batch/jet-job.sh 2760 $i 0.6 Fake -1 1;
done 

