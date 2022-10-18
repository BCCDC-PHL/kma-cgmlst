process kma_align {

    tag { sample_id }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_kma.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(read_1), path(read_2), val(scheme)

    output:
    tuple val(sample_id), path("${sample_id}_kma.csv")

    script:
    """
    ln -s ${scheme}.comp.b .
    ln -s ${scheme}.length.b .
    ln -s ${scheme}.name .
    ln -s ${scheme}.seq.b .
    
    kma \
      -t ${task.cpus} \
      -cge \
      -boot \
      -1t1 \
      -mem_mode \
      -and \
      -o ${sample_id}.kma \
      -t_db ${scheme} \
      -ipe ${read_1} ${read_2} \
      -tmp .

    head -n 1 ${sample_id}.kma.res | tr -d '#' | awk '{print tolower(\$0)}' > ${sample_id}_kma.csv
    tail -qn+2 ${sample_id}.kma.res | tr -d ' ' | tr \$'\\t' ',' >> ${sample_id}_kma.csv
    """
}
