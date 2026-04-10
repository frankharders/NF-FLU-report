process EXPORT_CONSENSUS {

    tag "$sampleName"

    publishDir "${params.outdir}/consensus_renamed", mode: 'copy'

    input:
    tuple val(sampleName), path(consensus_fasta)

    output:
    path("${sampleName}.renamed.consensus.fasta"), emit: fasta

    script:
    """
    python3 << 'EOF'
from pathlib import Path

sample = "${sampleName}"
fasta_in = Path("${consensus_fasta}")
fasta_out = Path(f"{sample}.renamed.consensus.fasta")

def normalize_fragment(header: str) -> str:
    h = header.strip().lstrip(">").split()[0]

    known = [
        "PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS",
        "M", "segment1", "segment2", "segment3", "segment4",
        "segment5", "segment6", "segment7", "segment8"
    ]

    for k in known:
        if k in h:
            if k == "segment1":
                return "PB2"
            elif k == "segment2":
                return "PB1"
            elif k == "segment3":
                return "PA"
            elif k == "segment4":
                return "HA"
            elif k == "segment5":
                return "NP"
            elif k == "segment6":
                return "NA"
            elif k == "segment7":
                return "MP"
            elif k == "segment8":
                return "NS"
            elif k == "M":
                return "MP"
            else:
                return k

    return h.replace(" ", "_")

with open(fasta_in) as fin:
    records = []
    current_header = None
    seq_lines = []

    for line in fin:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_header is not None:
                fragment = normalize_fragment(current_header)
                sequence = "".join(seq_lines).replace(" ", "")
                records.append((fragment, sequence))
            current_header = line
            seq_lines = []
        else:
            seq_lines.append(line)

    if current_header is not None:
        fragment = normalize_fragment(current_header)
        sequence = "".join(seq_lines).replace(" ", "")
        records.append((fragment, sequence))

with open(fasta_out, "w") as fout:
    for fragment, sequence in records:
        fout.write(f">{sample}_{fragment}\\n")
        fout.write(f"{sequence}\\n")

EOF
    """
}