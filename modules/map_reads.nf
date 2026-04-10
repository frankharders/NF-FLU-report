process MAP_READS {

    tag "$sampleName"

    conda "bioconda::minimap2 bioconda::samtools"

    input:
    tuple val(sampleName), val(platform), path(trimmed_reads)
    path consensus

    output:
    tuple val(sampleName), path("${sampleName}.sorted.bam"), path("${sampleName}.sorted.bam.bai"), emit: bam

    script:
    def preset = platform == "illumina" ? "sr" : "map-ont"
    """
    # Index consensus
    minimap2 -d ref.mmi ${consensus}

    # Map reads
    if [ "${platform}" == "illumina" ]; then
        minimap2 \
            -ax ${preset} \
            -t 4 \
            ref.mmi \
            ${trimmed_reads[0]} \
            ${trimmed_reads[1]} \
            | samtools sort -o ${sampleName}.sorted.bam
    else
        minimap2 \
            -ax ${preset} \
            -t 4 \
            ref.mmi \
            ${trimmed_reads[0]} \
            | samtools sort -o ${sampleName}.sorted.bam
    fi

    # Index bam
    samtools index ${sampleName}.sorted.bam

    # Rsync output
    mkdir -p ${params.outdir}/bam
    rsync -av ${sampleName}.sorted.bam     ${params.outdir}/bam/
    rsync -av ${sampleName}.sorted.bam.bai ${params.outdir}/bam/
    """
}
