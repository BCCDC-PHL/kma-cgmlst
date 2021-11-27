process parse_fastp_json {

    tag { sample_id }

    executor 'local'

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp_fastp.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(fastp_json)

    output:
    tuple val(sample_id), path("${sample_id}_fastp.csv")
    
    script:
    """
    parse_fastp_json.py \
      --fastp_json "${fastp_json}" \
      --sample_id "${sample_id}" \
      > ${sample_id}_fastp.csv
    """
}