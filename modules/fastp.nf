process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_fastp.{json,csv}"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastp.csv"), emit: csv

    script:
    """
    fastp \
      -t ${task.cpus} \
      --cut_tail \
      -i ${reads_1} \
      -I ${reads_2} \
      -o ${sample_id}_R1.trim.fastq.gz \
      -O ${sample_id}_R2.trim.fastq.gz \
      -j ${sample_id}_fastp.json

    parse_fastp_json.py \
      --fastp_json ${sample_id}_fastp.json \
      --sample_id "${sample_id}" \
      > ${sample_id}_fastp.csv
    """
}
