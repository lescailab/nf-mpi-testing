#!/usr/bin/env nextflow

/*

this small pipeline is based on the following implementation
of BWA based on MPI
https://github.com/bioinfo-pf-curie/mpiBWA/blob/master/docs/README.md#mpibwa

and attempts following the guidelines described here
https://www.nextflow.io/docs/latest/ignite.html#execution-with-mpi
*/
nextflow.enable.dsl=2

params.reads = '/home/lescai/tests/nextflow/small_data_human/*_R{1,2}_xxx.fastq.gz'
params.reference = '/home/lescai/tests/nextflow/short_refs/human_g1k_v37_decoy.small.fasta'
params.outdir = '/home/lescai/tests/nextflow/mpi/output'
params.publish_dir_mode = "copy"

process MPIBWA {

  tag "mpirun bwa"
  container "/home/lescai/singularity_folders/pull/mpiBWA.sif"

  publishDir "${params.outdir}/",
             mode: params.publish_dir_mode


  input:
  tuple val(sample), path(reads)
  path(reference)

  output:
  path("${sample}_alignment.sam")

  script:
  """
  mpiBWA mem \\
  -t ${task.cpus} \\
  -o "${sample}_alignment.sam" \\
  ${reference} \\
  $reads[0] $reads[1]
  """
}

workflow {

  inputReads = Channel.fromFilePairs(params.reads) 
  referenceFile = Channel.value(file(params.reference))

  MPIBWA(inputReads, referenceFile)

}