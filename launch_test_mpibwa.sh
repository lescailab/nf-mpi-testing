#!/bin/bash
#SBATCH --job-name=nfmpi
#SBATCH --output=/home/lescai/tests/nextflow/mpi/test_mpibwa_%j.log
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=16
#SBATCH --tasks-per-node=1

module load conda/anaconda3
module load singularity/3.3.0
module load openmpi/gcc/64/1.10.7

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/cm/shared/apps/anaconda/3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/cm/shared/apps/anaconda/3/etc/profile.d/conda.sh" ]; then
        . "/cm/shared/apps/anaconda/3/etc/profile.d/conda.sh"
    else
        export PATH="/cm/shared/apps/anaconda/3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate nextflow_21.04.0

## export NXF_MODE=ignite ### this doesn't make any difference in reality
export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
srun nextflow run \
-c /home/lescai/tests/nextflow/mpi/eos_mpi.conf \
-profile singularity \
/home/lescai/tests/nextflow/mpi/test_mpibwa.nf \
-process.executor ignite \ ### this makes all the difference - with this it works
-with-mpi