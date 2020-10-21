process combine_parsed_fastp_reports {
    executor 'local'
    publishDir "${params.outdir}/qc", pattern: "all.fastp_summary.tsv", mode: 'copy'

    input:
    path(parsed_fastp_reports)

    output:
    path("all.fastp_summary.tsv")
    
    script:
    first_parsed_report = parsed_fastp_reports[0]
    """
    head -n 1 ${first_parsed_report} > header.tsv
    tail -qn+2 ${parsed_fastp_reports} > data.tsv
    cat header.tsv data.tsv > all.fastp_summary.tsv
    """
}