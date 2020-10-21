process kma_result_to_mlst {
    tag { sample_id }
    executor 'local'
    publishDir "${params.outdir}", pattern: "${sample_id}_mlst.tsv", mode: 'copy'
    input:
    tuple val(sample_id), path(kma_result)

    output:
    tuple val(sample_id), path("${sample_id}_mlst.tsv")
    
    script:
    """
    kma_result_to_mlst.py \
      --res "${kma_result}" \
      --locus_allele_delimiter "_" \
      > ${sample_id}_mlst.tsv
    """
}