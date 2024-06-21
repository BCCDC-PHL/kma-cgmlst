#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { hash_files as hash_fastq_short } from './modules/hash_files.nf'
include { hash_files as hash_fastq_long }  from './modules/hash_files.nf'
include { fastp }                          from './modules/fastp.nf'
include { filtlong }                       from './modules/long_read_qc.nf'
include { nanoq as nanoq_pre_filter }      from './modules/long_read_qc.nf'
include { nanoq as nanoq_post_filter }     from './modules/long_read_qc.nf'
include { merge_nanoq_reports }            from './modules/long_read_qc.nf'
include { kma_align }                      from './modules/kma_align.nf'
include { kma_result_to_mlst }             from './modules/kma_result_to_mlst.nf'
include { count_called_alleles }           from './modules/count_called_alleles.nf'
include { pipeline_provenance }            from './modules/provenance.nf'
include { collect_provenance }             from './modules/provenance.nf'


if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

workflow {

    ch_workflow_metadata = Channel.value([
	workflow.sessionId,
	workflow.runName,
	workflow.manifest.name,
	workflow.manifest.version,
	workflow.start,
    ])

    if (params.samplesheet_input != 'NO_FILE') {
	ch_illumina_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['R1'], it['R2']]] }
	ch_nanopore_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['LONG']]] }.filter{ it -> it[1] != null }
    } else {
	ch_illumina_fastq = Channel.fromFilePairs( params.fastq_illumina_search_path, flat: true ).map{ it -> [it[0].split('_')[0], [it[1], it[2]]] }.unique{ it -> it[0] }
	ch_nanopore_fastq = Channel.fromPath( params.fastq_nanopore_search_path ).map{ it -> [it.getName().split('_')[0], [it]] }.unique{ it -> it[0] }
    }

    ch_scheme = Channel.fromPath( "${params.scheme}")
    
    main:
    ch_illumina_sample_ids = ch_illumina_fastq.map{ it -> it[0] }
    ch_nanopore_sample_ids = ch_nanopore_fastq.map{ it -> it[0] }
    ch_sample_ids = ch_illumina_sample_ids.concat(ch_nanopore_sample_ids).unique()

    hash_fastq_short(ch_illumina_fastq.combine(Channel.of("fastq-input")))
    hash_fastq_long(ch_nanopore_fastq.combine(Channel.of("fastq-input")))

    fastp(ch_illumina_fastq)

    nanoq_pre_filter(ch_nanopore_fastq.combine(Channel.of("pre_filter")))

    filtlong(ch_nanopore_fastq)

    nanoq_post_filter(filtlong.out.filtered_reads.combine(Channel.of("post_filter")))

    merge_nanoq_reports(nanoq_pre_filter.out.report.join(nanoq_post_filter.out.report))

    trimmed_reads = fastp.out.trimmed_reads.mix(filtlong.out.filtered_reads.map{ it -> [it[0], [it[1]]] })
    
    kma_align(trimmed_reads.combine(ch_scheme).view())

    kma_result_to_mlst(kma_align.out.res.combine(ch_scheme))

    count_called_alleles(kma_result_to_mlst.out.mlst)

    if (params.collect_outputs) {
	fastp.out.csv.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_fastp.csv", storeDir: params.outdir, keepHeader: true, sort: { it -> it.readLines()[1].split(',')[0] })

	count_called_alleles.out.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_called_allele_count.csv", storeDir: params.outdir, keepHeader: true, sort: { it -> it.readLines()[1].split(',')[0] })

	kma_result_to_mlst.out.mlst.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_cgmlst.csv", storeDir: params.outdir, keepHeader: true, sort: { it -> it.readLines()[1].split(',')[0] })
    }

    // Collect Provenance
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // At each step, we add another provenance file to the list using the << operator...
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_provenance = ch_sample_ids
    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_fastq_short.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(kma_align.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)
}
