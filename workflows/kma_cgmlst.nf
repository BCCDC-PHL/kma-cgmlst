nextflow.enable.dsl = 2

include { fastp } from '../modules/fastp.nf'
include { kma_align } from '../modules/kma_align.nf'
include { kma_result_to_mlst } from '../modules/kma_result_to_mlst.nf'

workflow kma_cgmlst {
    take:
      ch_fastq_input
      ch_scheme

    main:
      fastp (
        ch_fastq_input
      )

      kma_align(
        fastp.out[0]
	.combine(ch_scheme)
      )

      kma_result_to_mlst(
        kma_align.out
      )
}