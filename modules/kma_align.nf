process kma_align {

    tag { sample_id }
    publishDir "${params.outdir}", pattern: "${sample_id}.kma.res", mode: 'copy'
    input:
    tuple val(sample_id), path(read_1), path(read_2), val(scheme)

    output:
    tuple val(sample_id), path("${sample_id}.kma.res")

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
      -ipe ${read_1} ${read_2}
    """
}