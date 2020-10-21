nextflow.enable.dsl = 2

include { fastp } from '../modules/fastp.nf'
include { parse_fastp_json } from '../modules/parse_fastp_json.nf'
include { combine_parsed_fastp_reports } from '../modules/combine_parsed_fastp_reports.nf'
include { kma_align } from '../modules/kma_align.nf'
include { kma_result_to_mlst } from '../modules/kma_result_to_mlst.nf'

workflow kma_cgmlst {
    take:
      ch_fastq_input
      ch_scheme

    main:
      fastp(
        ch_fastq_input
      )

      parse_fastp_json(
        fastp.out[1]
      )

      combine_parsed_fastp_reports(
        parse_fastp_json.out.collect()
      )

      kma_align(
        fastp.out[0]
	.combine(ch_scheme)
      )

      kma_result_to_mlst(
        kma_align.out
      )
}