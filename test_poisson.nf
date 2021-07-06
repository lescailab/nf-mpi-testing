#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.outdir = '/home/lescai/tests/nextflow/mpi/poisson'
params.publish_dir_mode = "copy"

process POISSON {

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