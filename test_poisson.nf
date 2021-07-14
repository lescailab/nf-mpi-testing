#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.outdir = './poisson'
params.publish_dir_mode = "copy"
params.inputFile = '/home/lescai/tests/nextflow/mpi/demo_poisson_mod.py'

process POISSON {

  tag "fenics on ${pyfile}"
  container "/home/lescai/tests/nextflow/mpi/fenicsproject_2019.1.0.sif"

  publishDir "${params.outdir}",
             mode: params.publish_dir_mode

  input:
  path(pyfile)

  output:
  path("*vtu")
  path("*pvd")

  script:
  """
  python3 ${pyfile}
  """
}

workflow {

  inputFile = Channel.fromPath(params.inputFile) 

  POISSON(inputFile)

}