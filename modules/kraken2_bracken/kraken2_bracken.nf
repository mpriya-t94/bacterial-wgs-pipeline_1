process KRAKEN2_CLASSIFY {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/classification/kraken2", mode: 'copy'

    // Docker container for Kraken2
    container 'biocontainers/kraken2:2.1.3--pl5321h9f5acd7_0'

    input:
    tuple val(meta), path(reads)
    path kraken2_db from params.kraken2_db

    output:
    tuple val(meta.id), path("${meta.id}.kraken2.report"), emit: report
    tuple val(meta.id), path("${meta.id}.kraken2.mpa"), emit: mpa
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    def input_cmd = params.single_end ? "\"${reads}\"" : "--paired \"${reads[0]}\" \"${reads[1]}\""
    def report_zero_counts = params.report_zero_counts ? "--report-zero-counts" : ""
    """
    set -euo pipefail

    kraken2 \\
        --db "${params.kraken2_db}" \\
        --threads "${task.cpus}" \\
        "${report_zero_counts}" \\
        --report "${prefix}.kraken2.report" \\
        --confidence 0.1 \\
        "${input_cmd}" > /dev/null

    kreport2mpa.py \\
        -r "${prefix}.kraken2.report" \\
        -o "${prefix}.kraken2.mpa" \\
        --display-header

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}

process BRACKEN_ABUNDANCE {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/classification/bracken", mode: 'copy'

    // Docker container for Bracken
    container 'biocontainers/bracken:2.8--py39hc16433a_0'

    input:
    tuple val(meta.id), path(kraken2_report) from KRAKEN2_CLASSIFY.out.report
    path bracken_db from params.bracken_db

    output:
    tuple val(meta.id), path("${meta.id}.bracken"), emit: bracken
    tuple val(meta.id), path("${meta.id}.bracken.mpa"), emit: bracken_mpa
    tuple val(meta.id), path("${meta.id}.bracken.report"), emit: bracken_report
    path "versions.yml", emit: versions
    
    when:
    !params.skip_bracken

    script:
    def prefix = "${meta.id}"
    def read_length = params.read_length ?: 150
    def level = params.bracken_level ?: 'S'
    """
    set -euo pipefail
    bracken \\
        -d "${bracken_db}" \\
        -i "${kraken2_report}" \\
        -o "${meta.id}.${level}.bracken" \\
        -w "${meta.id}.${level}.bracken.report" \\
        -r "${read_length}" \\
        -l "${level}" \\
        -t 10

    kreport2mpa.py \\
        -r "${meta.id}.${level}.bracken" \\
        -o "${meta.id}.bracken.mpa" \\
        --display-header

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: \$(echo \$(bracken --version 2>&1) | sed 's/^.*Bracken version //; s/ .*\$//')
    END_VERSIONS
    """
}

process COMBINE_KRAKEN2_REPORTS {
    label "process_single"
    publishDir "${params.outdir}/combined/kraken2", mode: 'copy'

    container 'biocontainers/kraken2:2.1.3--pl5321h9f5acd7_0'

    input:
    path kraken2_reports from KRAKEN2_CLASSIFY.out.report.collect()

    output:
    path "combined_kraken2.report", emit: combined_report

    when:
    !params.skip_combined_reports

    script:
    """
    set -euo pipefail
    combine_kraken_reports.py \\
        -r "${kraken2_reports.join(' ')}" \\
        -o combined_kraken2.report
    """
}

process COMBINE_KRAKEN2_MPA {
    label "process_single"
    publishDir "${params.outdir}", mode: "copy"

    container 'biocontainers/kraken2:2.1.3--pl5321h9f5acd7_0'

    input:
    path(kraken2_mpa).collect() from KRAKEN2_CLASSIFY.out.mpa

    output:
    path "combined.kraken2.mpa", emit: combined_mpa

    when:
    !params.skip_combining

    script:
    """
    set -euo pipefail
    combine_mpa.py \\
        -i "${kraken2_mpa.join(' ')}" \\
        -o combined.kraken2.mpa
    """
}

process COMBINE_BRACKEN_REPORTS {
    label "process_single"
    publishDir "${params.outdir}/combined/bracken", mode: 'copy'

    container 'biocontainers/bracken:2.8--py39hc16433a_0'

    input:
    path(bracken_reports).collect() from BRACKEN_ABUNDANCE.out.bracken_report.collect()

    output:
    path "combined_bracken.report", emit: combined_bracken_report

    when:
    !params.skip_bracken && !params.skip_combined_reports

    script:
    """
    set -euo pipefail
    combine_bracken_reports.py \\
        -r "${bracken_reports.join(' ')}" \\
        -o combined_bracken.report
    """
}

process COMBINE_BRACKEN_MPA {
    tag "${level}"
    label 'process_single'
    publishDir "${params.outdir}/combined_bracken_mpa", mode: 'copy'
    
    container 'biocontainers/kraken2:2.1.3--pl5321h9f5acd7_0'
    
    input:
    tuple val(level), path(bracken_mpa).collect() from BRACKEN_ABUNDANCE.out.bracken_mpa.collect()
    
    output:
    path "combined.${level}.bracken.mpa", emit: combined_bracken_mpa
    
    when:
    !params.skip_combining && !params.skip_bracken
    
    script:
    """
    set -euo pipefail
    combine_mpa.py \\
        -i "${bracken_mpa.join(' ')}" \\
        -o combined.${level}.bracken.mpa
    """
}
