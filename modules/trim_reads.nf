process TRIM_READS {

    tag "$sampleName"

    conda "bioconda::fastp=0.23.4"

    input:
    tuple val(sampleName), val(platform), path(reads)

    output:
    tuple val(sampleName), val(platform), path("${sampleName}_trimmed*.fastq.gz"), emit: trimmed_reads
    path "${sampleName}_fastp.json", emit: fastp_json
    path "${sampleName}_fastp.html", emit: fastp_html

    script:
    if (platform == "illumina")
        """
        fastp \
            --in1 ${reads[0]} \
            --in2 ${reads[1]} \
            --out1 ${sampleName}_trimmed_R1.fastq.gz \
            --out2 ${sampleName}_trimmed_R2.fastq.gz \
            --json ${sampleName}_fastp.json \
            --html ${sampleName}_fastp.html \
            --thread 4 \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 50

        mkdir -p ${params.outdir}/trimmed
        rsync -av ${sampleName}_trimmed_R1.fastq.gz ${params.outdir}/trimmed/
        rsync -av ${sampleName}_trimmed_R2.fastq.gz ${params.outdir}/trimmed/
        rsync -av ${sampleName}_fastp.json          ${params.outdir}/trimmed/
        rsync -av ${sampleName}_fastp.html          ${params.outdir}/trimmed/
        """
    else if (platform == "nanopore")
        """
        chopper \
            --input ${reads[0]} \
            --quality 8 \
            --minlength 200 \
            | gzip > ${sampleName}_trimmed.fastq.gz

        echo '{}' > ${sampleName}_fastp.json
        echo '<html></html>' > ${sampleName}_fastp.html

        mkdir -p ${params.outdir}/trimmed
        rsync -av ${sampleName}_trimmed.fastq.gz ${params.outdir}/trimmed/
        rsync -av ${sampleName}_fastp.json       ${params.outdir}/trimmed/
        rsync -av ${sampleName}_fastp.html       ${params.outdir}/trimmed/
        """
}
