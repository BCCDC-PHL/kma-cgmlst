process kma_result_to_mlst {
    tag { sample_id }

    executor 'local'

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_{cgmlst,locus_qc}.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(kma_result), path(schemeName)

    output:
    tuple val(sample_id), path("${sample_id}_cgmlst.csv"), emit: mlst
    tuple val(sample_id), path("${sample_id}_locus_qc.csv"), emit: mlst_qc
    
    script:
    """
   
    kma_result_to_mlst.py \
      "${kma_result}" \
      --alleles "${schemeName}" \
      --sample-id "${sample_id}" \
      --locus-allele-delimiter "_" \
      --min-identity ${params.min_identity} \
      --min-coverage ${params.min_coverage} \
      -o ${sample_id}_cgmlst.csv \
      > ${sample_id}_locus_qc.csv
    """
}
