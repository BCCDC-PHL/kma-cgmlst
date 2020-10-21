process fastqc {
    tag { sample_id }
    label 'cpu4'

    input:
    tuple val(grouping_key), path(fastq)

    output:
    tuple val(sample_id), path("${sample_id}.kma.res")

    script:
    if (grouping_key =~ '_S[0-9]+_') {
      sample_id = grouping_key.split("_S[0-9]+_")[0]
    } else {
      sample_id = grouping_key.split("_")[0]
    }

    """
    fastqc \
      -t ${task.cpus} \
      -f fastq \
      ${fastq}
    unzip 
    """
}