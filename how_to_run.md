# How to run test on EOS

## Have your environment ready

You can get detailed information about conda environments on EOS here, and how to run Nextflow on EOS [here](https://lescailab.unipv.it/guides/eos_guide/use_nextflow.html).

In short:

First of all, load EOS modules you will need

```
module load conda/anaconda3
```

you might have to intialise conda in the very same bash shell where you are by running

```
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
```

Then you need to create a conda environment where you can install nextflow

```
conda create --name nextflow_21.04.0 -c bioconda nextflow=21.04.0
```

Now you're ready.

## Request resources for nextflow process

You will run the nextflow master from an interactive job, by launching

```
srun -c 1 --mem 5g --job-name nf-master --pty bash
```

Once you get resources assigned, you need to activate all your tools

```
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
```

You are now ready to run.

Make sure you have copied the following files in the folder you are located in:

- eos.config
- test_poisson_mpi.nf
- demo_poisson_mod.py

Now you can launch the test workflow by typing:

```
nextflow run \
-c eos.config \
-profile eos \
test_poisson_mpi.nf \
--inputFile ./demo_poisson_mod.py
```

When the analysis is finished you will find a new folder in the directory you're located in named ```poisson``.
Inside that folder, the expected results of the run.