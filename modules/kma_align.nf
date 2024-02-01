process kma_align {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_kma*.{c,t}sv", mode: 'copy'

    input:
    tuple val(sample_id), path(read_1), path(read_2), val(scheme)

    output:
    tuple val(sample_id), path("${sample_id}_kma.csv"), emit: res
    tuple val(sample_id), path("${sample_id}_kma_mapstat.tsv"), emit: mapstat

    script:
    """
    ln -s ${scheme}.comp.b .
    ln -s ${scheme}.length.b .
    ln -s ${scheme}.name .
    ln -s ${scheme}.seq.b .
    
    kma \
      -t ${task.cpus} \
      -ef \
      -cge \
      -boot \
      -1t1 \
      -mem_mode \
      -and \
      -o ${sample_id}.kma \
      -t_db ${scheme} \
      -ipe ${read_1} ${read_2} \
      -tmp .

    head -n 1 ${sample_id}.kma.res | tr -d '#' | awk '{print tolower(\$0)}' | tr \$'\\t' ',' > ${sample_id}_kma.csv
    tail -qn+2 ${sample_id}.kma.res | tr -d ' ' | tr \$'\\t' ',' >> ${sample_id}_kma.csv

    mv ${sample_id}.kma.mapstat ${sample_id}_kma_mapstat.tsv
    """
}
