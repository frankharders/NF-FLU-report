process COVERAGE_PLOT {

    tag "$sampleName"

    conda "conda-forge::matplotlib conda-forge::pandas bioconda::pysam"

    input:
    tuple val(sampleName), path(vcf), path(tbi)
    tuple val(sampleName2), path(depth)
    path gff
    path mutation_check
    path plot_script

    output:
    path "${sampleName}_coverage.pdf", emit: pdf

    script:
    """
    python3 ${plot_script} \
        ${sampleName} \
        ${vcf} \
        ${depth} \
        ${gff} \
        ${mutation_check}

    mkdir -p ${params.outdir}/plots
    rsync -av ${sampleName}_coverage.pdf ${params.outdir}/plots/
    """
}
