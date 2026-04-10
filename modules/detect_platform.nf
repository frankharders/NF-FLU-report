process DETECT_PLATFORM {

    input:
    path samplesheet

    output:
    env PLATFORM, emit: platform

    script:
    """
    HEADER=\$(head -1 ${samplesheet})

    if echo "\$HEADER" | grep -q "fastq_2"; then
        echo "illumina" > platform.txt
    elif echo "\$HEADER" | grep -q "reads"; then
        echo "nanopore" > platform.txt
    else
        echo "unknown" > platform.txt
    fi

    PLATFORM=\$(cat platform.txt)
    echo "Detected platform: \$PLATFORM"
    """
}
