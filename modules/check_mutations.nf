process CHECK_MUTATIONS {

    tag "$sampleName"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sampleName), path(faa), path(flumut_parsed)

    output:
    path "${sampleName}_mutation_check.tsv", emit: mutation_check

    script:
    """
    python3 << 'EOF'
import re
import csv

# --- Parse FAA: gen -> aminozuursequentie ---
faa_seqs = {}
current_gene = None
current_seq  = []

with open("${faa}", "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_gene and current_seq:
                faa_seqs[current_gene] = "".join(current_seq)
            parts     = line.split("|")
            feat_type = parts[1] if len(parts) > 1 else ""
            gene_name = parts[2] if len(parts) > 2 else ""
            current_gene = gene_name if feat_type == "CDS" else None
            current_seq  = []
        else:
            if current_gene:
                current_seq.append(line)

    if current_gene and current_seq:
        faa_seqs[current_gene] = "".join(current_seq)

# --- Gen mapping ---
def resolve_gene(segment_field):
    mapping = {
        "HA1-5" : "HA",
        "HA2-5" : "HA",
        "NS-1"  : "NS1",
        "NS1"   : "NS1",
        "M1"    : "M1",
        "M2"    : "M2",
        "NP"    : "NP",
        "NA"    : "NA",
        "PB2"   : "PB2",
        "PB1"   : "PB1",
        "PA"    : "PA",
    }
    return mapping.get(segment_field, segment_field)

# --- Verwerk flumut_parsed.tsv ---
header = ["Sample", "Segment", "Gene", "Position", "WT_aa", "Mut_aa",
          "Observed_aa", "Mutation_Present", "Effects", "Subtype", "Literature"]

with open("${sampleName}_mutation_check.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(header)

    with open("${flumut_parsed}", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["Sample"].strip()
            if sample != "${sampleName}":
                continue

            segment  = row["Segment"].strip()
            wt_aa    = row["WT_aa"].strip()
            position = row["Position"].strip()
            mut_aa   = row["Mut_aa"].strip()
            effects  = row["Effects"].strip()
            subtype  = row["Subtype"].strip()
            literature = row["Literature"].strip()

            gene = resolve_gene(segment)
            pos  = int(position)
            seq  = faa_seqs.get(gene, None)

            if seq is None:
                observed = "gene_not_found"
                present  = "UNKNOWN"
            elif pos < 1 or pos > len(seq):
                observed = "pos_out_of_range"
                present  = "UNKNOWN"
            else:
                observed = seq[pos - 1]
                present  = "YES" if observed == mut_aa else "NO"

            writer.writerow([
                sample, segment, gene, position,
                wt_aa, mut_aa, observed,
                present, effects, subtype, literature
            ])

EOF
    """
}
