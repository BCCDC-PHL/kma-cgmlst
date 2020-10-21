process parse_fastp_json {
    tag { sample_id }
    executor 'local'
    input:
    tuple val(sample_id), path(fastp_json)

    output:
    tuple val(sample_id), path("${sample_id}_fastp.tsv")
    
    script:
    """
    parse_fastp_json.py \
      --fastp_json "${fastp_json}" \
      > ${sample_id}_fastp.tsv
    """
}