import nextflow.util.BlankSeparatedList
import nextflow.processor.TaskPath

process kma_align {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_kma*.{c,t}sv", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(scheme)

    output:
    tuple val(sample_id), path("${sample_id}_kma.csv"), emit: res
    tuple val(sample_id), path("${sample_id}_kma_mapstat.tsv"), emit: mapstat
    tuple val(sample_id), path("${sample_id}_kma_align_provenance.yml"), emit: provenance

    script:
    reads_input = (reads instanceof nextflow.processor.TaskPath) ? "-i ${reads[0]}" : "-ipe ${reads[0]} ${reads[1]}"
    """
    printf -- "- process_name: kma_align\\n"       >> ${sample_id}_kma_align_provenance.yml
    printf -- "  tools:\\n"                        >> ${sample_id}_kma_align_provenance.yml
    printf -- "    - tool_name: kma\\n"            >> ${sample_id}_kma_align_provenance.yml
    printf -- "      tool_version: \$(kma -v 2>&1 | cut -d '-' -f 2)\\n" >> ${sample_id}_kma_align_provenance.yml
    printf -- "      parameters:\\n"               >> ${sample_id}_kma_align_provenance.yml
    printf -- "        - parameter: -ef\\n"        >> ${sample_id}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${sample_id}_kma_align_provenance.yml
    printf -- "        - parameter: -cge\\n"       >> ${sample_id}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${sample_id}_kma_align_provenance.yml
    printf -- "        - parameter: -boot\\n"      >> ${sample_id}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${sample_id}_kma_align_provenance.yml
    printf -- "        - parameter: -1t1\\n"       >> ${sample_id}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${sample_id}_kma_align_provenance.yml
    printf -- "        - parameter: -mem_mode\\n"  >> ${sample_id}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${sample_id}_kma_align_provenance.yml
    printf -- "        - parameter: -t_db\\n"      >> ${sample_id}_kma_align_provenance.yml
    printf -- "          value: ${scheme}\\n"      >> ${sample_id}_kma_align_provenance.yml
    printf -- "        - parameter: -and\\n"       >> ${sample_id}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${sample_id}_kma_align_provenance.yml

    # ln -s ${scheme}.comp.b .
    # ln -s ${scheme}.length.b .
    # ln -s ${scheme}.name .
    # ln -s ${scheme}.seq.b .
    
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
	${reads_input} \
	-tmp .

    head -n 1 ${sample_id}.kma.res | tr -d '#' | awk '{print tolower(\$0)}' | tr \$'\\t' ',' > ${sample_id}_kma.csv
    tail -qn+2 ${sample_id}.kma.res | tr -d ' ' | tr \$'\\t' ',' >> ${sample_id}_kma.csv

    mv ${sample_id}.kma.mapstat ${sample_id}_kma_mapstat.tsv
    """
}
