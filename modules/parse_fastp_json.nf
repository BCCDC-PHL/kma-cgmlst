process parse_fastp_json {
    tag { sample_id }
    executor 'local'
    input:
    tuple val(sample_id), path(fastp_json)

    output:
    path("${sample_id}_fastp.tsv")
    
    script:
    """
    parse_fastp_json.py \
      --fastp_json "${fastp_json}" \
      --sample_id "${sample_id}" \
      > ${sample_id}_fastp.tsv
    """
}