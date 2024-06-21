process TAG_BARCODE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/sccore.yml"
    container "singleron-rd/sccore:1.0.0"

    input:
    tuple val(meta), path(reads), path(match_barcode)
    path tag_barcode_fasta
    val r2_pattern

    output:
    tuple val(meta), path("*.json"),  emit: json
    tuple val(meta), path("*.csv*"),  emit: csv
    path  "versions.yml" , emit: versions

    script:
    // separate forward from reverse pairs
    """
    tag_barcode.py \\
        --sample ${meta.id} \\
        --fq ${reads} \\
        --tag_barcode_fasta ${tag_barcode_fasta} \\
        --match_barcode ${match_barcode} \\
        --r2_pattern ${r2_pattern}
  
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccore: \$(format-bc --version)
    END_VERSIONS
    """
}