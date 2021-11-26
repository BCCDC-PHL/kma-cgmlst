#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp } from './modules/fastp.nf'
include { parse_fastp_json } from './modules/parse_fastp_json.nf'
include { combine_parsed_fastp_reports } from './modules/combine_parsed_fastp_reports.nf'
include { kma_align } from './modules/kma_align.nf'
include { kma_result_to_mlst } from './modules/kma_result_to_mlst.nf'
include { count_called_alleles } from './modules/count_called_alleles.nf'

if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

workflow {
    Channel.fromFilePairs( "${params.fastq_input}/*_{1,2}*.fastq.gz", type: 'file', maxDepth: 1 ).set{ ch_fastq_input }
    Channel.fromPath( "${params.scheme}").set{ ch_scheme }
    
    main:
      fastp(ch_fastq_input)

      parse_fastp_json(fastp.out[1])

      combine_parsed_fastp_reports(parse_fastp_json.out.collect())

      kma_align(fastp.out[0].combine(ch_scheme))

      kma_result_to_mlst(kma_align.out.combine(ch_scheme))

      count_called_alleles(kma_result_to_mlst.out)
}