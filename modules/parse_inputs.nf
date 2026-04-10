process PARSE_INPUTS {

    tag "$sampleName"

    input:
    tuple val(sampleName), path(cleavage), path(genin2)

    output:
    path "parsed_${sampleName}.txt", emit: report_input

    script:
    """
    LINE=\$(grep "^${sampleName}" ${genin2})
    GENOTYPE=\$(echo "\$LINE"    | awk -F'\t' '{print \$2}')
    SUBGENOTYPE=\$(echo "\$LINE" | awk -F'\t' '{print \$3}')

    CLINE=\$(grep "^${sampleName}" ${cleavage})
    CLEAVAGE_SEQ=\$(echo "\$CLINE" | awk -F'\t' '{print \$2}')
    RESIDUE_TYPE=\$(echo "\$CLINE" | awk -F'\t' '{print \$3}')
    PATHOGENICITY=\$(echo "\$CLINE" | awk -F'\t' '{print \$5}')

    echo -e "${sampleName}\t\$GENOTYPE\t\$SUBGENOTYPE\t\$CLEAVAGE_SEQ\t\$RESIDUE_TYPE\t\$PATHOGENICITY" \
        > parsed_${sampleName}.txt
    """
}
