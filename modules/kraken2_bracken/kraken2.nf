process KRAKEN2_CLASSIFY {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/classification/kraken2", mode: 'copy'

    container 'staphb/kraken2:2.1.2'

    input:
    tuple val(meta), path(reads)
    path kraken2_db

    output:
    tuple val(meta), path("${meta.id}.kraken2.report"), emit: report
    tuple val(meta), path("${meta.id}.kraken2.output"), emit: output
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    def input_cmd = (reads instanceof List && reads.size() == 2) ? 
        "--paired ${reads[0]} ${reads[1]}" : 
        "${reads}"
    
    """
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --report ${prefix}.kraken2.report \\
        --output ${prefix}.kraken2.output \\
        --confidence 0.1 \\
        ${input_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}

process PARSE_KRAKEN2 {
    tag "$meta.id"
    publishDir "${params.outdir}/${meta.id}/classification/kraken2", mode: 'copy'

    container 'python:3.9-slim'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("${meta.id}_top_organism.txt"), emit: top_organism
    tuple val(meta), env(ORGANISM), emit: organism_meta

    script:
    def prefix = "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import re
    
    def parse_kraken_report(report_file):
        """Extract top genus-level classification from Kraken2 report"""
        top_organism = "Unknown"
        max_percentage = 0.0
        
        with open(report_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\\t')
                if len(fields) >= 6:
                    percentage = float(fields[0])
                    rank_code = fields[3].strip()
                    tax_name = fields[5].strip()
                    
                    # Look for genus level (G) or species level (S)
                    if rank_code in ['G', 'S'] and percentage > max_percentage:
                        max_percentage = percentage
                        top_organism = tax_name
        
        return top_organism, max_percentage
    
    organism, percentage = parse_kraken_report("${report}")
    
    # Write to file
    with open("${prefix}_top_organism.txt", 'w') as out:
        out.write(f"Sample: ${meta.id}\\n")
        out.write(f"Top Organism: {organism}\\n")
        out.write(f"Percentage: {percentage}%\\n")
    
    # Export for channel
    print(f"ORGANISM={organism}")
    """
}
