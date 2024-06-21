process EXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/sccore.yml"
    container "singleron-rd/sccore:1.0.0"

    input:
    tuple val(meta), path(reads, stageAs: "?/*")
    path assets_dir
    val protocol

    output:
    tuple val(meta), path("${meta.id}_R2.fq.gz"),  emit: reads
    tuple val(meta), path("*.json"),  emit: json
    path  "versions.yml" , emit: versions

    script:
    // separate forward from reverse pairs
    def (r1,r2) = reads.collate(2).transpose()
    """
    extract.py \\
        --sample ${meta.id} \\
        --fq1 ${r1.join( "," )} \\
        --fq2 ${r2.join( "," )} \\
        --assets_dir ${assets_dir} \\
        --protocol ${protocol} 
   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccore: \$(format-bc --version)
    END_VERSIONS
    """
}