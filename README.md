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

To be continued...