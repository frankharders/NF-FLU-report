process CALL_VARIANTS {

    tag "$sampleName"

    conda "bioconda::bcftools bioconda::samtools"

    input:
    tuple val(sampleName), path(bam), path(bai)
    path consensus

    output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.tbi"), emit: vcf
    tuple val(sampleName), path("${sampleName}.depth.tsv"), emit: depth

    script:
    """
    bcftools mpileup \
        --fasta-ref ${consensus} \
        --annotate FORMAT/AD,FORMAT/DP \
        --max-depth 10000 \
        --min-BQ 20 \
        ${bam} \
    | bcftools call \
        --multiallelic-caller \
        --output-type z \
        --output ${sampleName}.raw.vcf.gz

    bcftools index --tbi ${sampleName}.raw.vcf.gz

    bcftools +fill-tags ${sampleName}.raw.vcf.gz -- -t AF \
    | bcftools view \
        --include 'AF>=0.01' \
        --output-type z \
        --output ${sampleName}.vcf.gz

    bcftools index --tbi ${sampleName}.vcf.gz

    samtools depth \
        -a \
        -d 10000 \
        ${bam} \
        > ${sampleName}.depth.tsv

    mkdir -p ${params.outdir}/variants
    rsync -av ${sampleName}.vcf.gz     ${params.outdir}/variants/
    rsync -av ${sampleName}.vcf.gz.tbi ${params.outdir}/variants/
    rsync -av ${sampleName}.depth.tsv  ${params.outdir}/variants/
    """
}
