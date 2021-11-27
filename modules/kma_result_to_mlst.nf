process kma_result_to_mlst {
    tag { sample_id }

    executor 'local'

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_{cgmlst,locus_qc}.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(kma_result), val(scheme)

    output:
    tuple val(sample_id), path("${sample_id}_cgmlst.csv"), path("${sample_id}_locus_qc.csv")
    
    script:
    """
    ln -s ${scheme}.name .
    kma_result_to_mlst.py \
      "${kma_result}" \
      --alleles ${scheme}.name \
      --sample-id "${sample_id}" \
      --locus-allele-delimiter "_" \
      -o ${sample_id}_cgmlst.csv \
      > ${sample_id}_locus_qc.csv
    """
}