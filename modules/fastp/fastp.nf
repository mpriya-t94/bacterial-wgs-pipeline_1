process FASTP {
    tag "$meta.id"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    // Docker container for FASTP
    container 'quay.io/biocontainers/fastp:1.0.1--heae3180_0'

    input:
    tuple val(meta), path(reads), path(fastqc_zip)

    output:
    tuple val(meta), path("*_trimmed_R*.fastq.gz"), emit: trimmed_reads
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.json"), emit: json

    script:
    def prefix = "${meta.id}"
    def r1 = reads[0]
    def r2 = reads[1]

    """
    unzip -q ${fastqc_zip}

    QUALITY_LINE=\$(unzip -p ${fastqc_zip} */summary.txt | grep "Per base sequence quality")

    QUALITY_THRESHOLD=20
    if echo "\$QUALITY_LINE" | grep -q "FAIL"; then
        QUALITY_THRESHOLD=25
    elif echo "\$QUALITY_LINE" | grep -q "WARN"; then
        QUALITY_THRESHOLD=22
    fi

    fastp \\
        --in1 ${r1} \\
        --in2 ${r2} \\
        --out1 ${prefix}_trimmed_R1.fastq.gz \\
        --out2 ${prefix}_trimmed_R2.fastq.gz \\
        --html ${prefix}_fastp.html \\
        --json ${prefix}_fastp.json \\
        --thread ${task.cpus} \\
        --qualified_quality_phred \$QUALITY_THRESHOLD \\
        --detect_adapter_for_pe \\
        --length_required 50
    """

}