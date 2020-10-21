#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { kma_cgmlst } from './workflows/kma_cgmlst.nf'

if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

workflow {
    Channel.fromFilePairs( "${params.fastq_input}/*_R{1,2}*.fastq.gz", type: 'file', maxDepth: 1 ).set{ ch_fastq }
    Channel.fromPath( "${params.scheme}").set{ ch_scheme }
    
    main:
      kma_cgmlst(
	ch_fastq,
	ch_scheme,
      )
}