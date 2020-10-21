process fastp {
    tag { sample_id }
    label 'cpu4'

    input:
    tuple val(grouping_key), path(fastq)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz")
    tuple val(sample_id), path("${sample_id}.fastp.json")
    

    script:
    if (grouping_key =~ '_S[0-9]+_') {
      sample_id = grouping_key.split("_S[0-9]+_")[0]
    } else {
      sample_id = grouping_key.split("_")[0]
    }
    read_1 = fastq[0]
    read_2 = fastq[1]
    """
    fastp \
      -t ${task.cpus} \
      -i ${read_1} \
      -I ${read_2} \
      -o ${sample_id}_R1.trim.fastq.gz \
      -O ${sample_id}_R2.trim.fastq.gz \
      -j ${sample_id}.fastp.json
    """
}