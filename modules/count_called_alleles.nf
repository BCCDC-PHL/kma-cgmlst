process count_called_alleles {

    tag { sample_id }

    executor 'local'

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_called_allele_count.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(cgmlst)

    output:
    tuple val(sample_id), path("${sample_id}_called_allele_count.csv")
    
    script:
    """
    count_called_alleles.py ${cgmlst} > ${sample_id}_called_allele_count.csv
    """
}