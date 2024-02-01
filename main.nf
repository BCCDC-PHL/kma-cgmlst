#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp }                        from './modules/fastp.nf'
include { kma_align }                    from './modules/kma_align.nf'
include { kma_result_to_mlst }           from './modules/kma_result_to_mlst.nf'
include { count_called_alleles }         from './modules/count_called_alleles.nf'

if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

workflow {
    if (params.samplesheet_input != 'NO_FILE') {
	ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true)
    } else {
	ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    }

    ch_scheme = Channel.fromPath( "${params.scheme}")
    
    main:
    fastp(ch_fastq)

    kma_align(fastp.out.trimmed_reads.combine(ch_scheme))

    kma_result_to_mlst(kma_align.out.res.combine(ch_scheme))

    count_called_alleles(kma_result_to_mlst.out.mlst)

    if (params.collect_outputs) {
	fastp.out.csv.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_fastp.csv", storeDir: params.outdir, keepHeader: true, sort: { it -> it.readLines()[1].split(',')[0] })

	count_called_alleles.out.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_called_allele_count.csv", storeDir: params.outdir, keepHeader: true, sort: { it -> it.readLines()[1].split(',')[0] })

	kma_result_to_mlst.out.mlst.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_cgmlst.csv", storeDir: params.outdir, keepHeader: true, sort: { it -> it.readLines()[1].split(',')[0] })
    }
}
