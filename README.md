# Nextflow MPI Tests

This is a test repository, with code used to assess use-cases and performance of Nextflow with openMPI tools.


## Key elements


As described in [here](https://github.com/nextflow-io/nextflow/issues/118), Nextflow doesn't package *ignite* dependencies, and therefore it needs to be informed at launch of the need to download it.
This important information should be added to the [documentation](https://www.nextflow.io/docs/latest/ignite.html#linux-slurm-launcher), but it is worth highlighting here.

The launch script has to include the option ```-process.executor ignite```, as in:

```
export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
srun nextflow run \
-c eos_mpi.conf \
-profile singularity \
test_mpibwa.nf \
-process.executor ignite \
-with-mpi
```

## Use cases


### FenicsProject

A solution is described [here](https://bitbucket.org/fenics-project/docker/src/master/dockerfiles/build-singularity-images.sh) to build a Singularity image.

However, we have tested first a conversion from Docker to Singularity in an environment with sudo rights, simply by:

```
singularity pull docker://docker pull quay.io/fenicsproject:stable
```

Then we test if it works, by using one of the demo scripts packed inside the image

```
singularity exec fenicsproject_2019.1.0.sif cp /usr/local/share/dolfin/demo/python/documented/poisson/demo_poisson.py .
```

We have modified the script in order to remove connection to display.
Then we execute:

```
singularity exec fenicsproject_2019.1.0.sif mpirun -np 12 python3 demo_poisson_mod.py
```

And it seems to work by broadcasting 12 processes, writing 12 .vtu files, one .pvtu and one output file ```poisson.pvd```.



To be continued...