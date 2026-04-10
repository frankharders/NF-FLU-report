process MAKE_REPORT {

    publishDir params.outdir, mode: 'copy'

    input:
    path parsed_files
    path mutation_files

    output:
    path "final_report.tsv"

    script:
    """
    python3 << 'EOF'
import csv
import glob

# --- Laad genotype/cleavage info per sample ---
sample_info = {}

for f in glob.glob("parsed_*.txt"):
    with open(f, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                continue
            sample        = cols[0]
            genotype      = cols[1]
            subgenotype   = cols[2]
            cleavage_seq  = cols[3]
            residue_type  = cols[4]
            pathogenicity = cols[5]
            sample_info[sample] = [genotype, subgenotype, cleavage_seq, residue_type, pathogenicity]

# --- Schrijf eindrapport ---
header = [
    "Sample", "Genotype", "Sub-genotype",
    "Cleavage Sequence", "Residue Type", "Pathogenicity",
    "Segment", "Gene", "Position", "WT_aa", "Mut_aa",
    "Observed_aa", "Mutation_Present", "Effects", "Subtype", "Literature"
]

with open("final_report.tsv", "w", newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(header)

    for f in sorted(glob.glob("*_mutation_check.tsv")):
        with open(f, "r") as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                sample = row['Sample']
                info   = sample_info.get(sample, ["NA"] * 5)
                writer.writerow([
                    sample,
                    info[0], info[1], info[2], info[3], info[4],
                    row['Segment'],
                    row['Gene'],
                    row['Position'],
                    row['WT_aa'],
                    row['Mut_aa'],
                    row['Observed_aa'],
                    row['Mutation_Present'],
                    row['Effects'],
                    row['Subtype'],
                    row['Literature']
                ])

EOF
    """
}
