process CHECK_CLEAVAGESITE {

    tag "$sampleName"

    conda "conda-forge::pandas bioconda::pysam"

    publishDir "${params.outdir}/cleavage_site", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(bai), path(gff)

    output:
    path("${sampleName}.cleavage_haplotypes.tsv"),         emit: haplotypes
    path("${sampleName}.cleavage_haplotypes.all.tsv"),     emit: haplotypes_all
    path("${sampleName}.cleavage_haplotypes.summary.tsv"), emit: summary

    script:
    def cs_min_mapq = params.containsKey('cs_min_mapq') ? params.cs_min_mapq : 20
    def cs_min_freq = params.containsKey('cs_min_freq') ? params.cs_min_freq : 0.01

    """
    python3 << 'EOF'
import pysam
import pandas as pd

sample = "${sampleName}"
bam_file = "${bam}"
gff_file = "${gff}"

MIN_MAPQ = int(${cs_min_mapq})
MIN_FREQ = float(${cs_min_freq})

# -----------------------------
# Parse GFF: HA cleavage site
# -----------------------------
cs_start = None
cs_end   = None
contig   = None

with open(gff_file) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\\n").split("\\t")
        if len(cols) < 9:
            continue

        seqid, source, feat, start, end, score, strand, phase, attrs = cols

        if feat == "misc_feature" and "gene=HA" in attrs and "cleavage site" in attrs:
            cs_start = int(start)
            cs_end   = int(end)
            contig   = seqid
            break

if cs_start is None or cs_end is None or contig is None:
    raise Exception(f"Cleavage site for HA not found in GFF: {gff_file}")

# -----------------------------
# Codon table
# -----------------------------
codon_table = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def translate(nt_seq):
    seq = nt_seq.replace("-", "").upper()
    aa = []
    for i in range(0, len(seq) - 2, 3):
        aa.append(codon_table.get(seq[i:i+3], "X"))
    return "".join(aa)

def full_span(read, start_1based, end_1based):
    return read.reference_start <= (start_1based - 1) and read.reference_end >= end_1based

bam = pysam.AlignmentFile(bam_file, "rb")

motifs = []
reads_seen_window = 0
reads_pass_mapq = 0
reads_pass_fullspan = 0
reads_with_deletion = 0

for read in bam.fetch(contig, cs_start - 1, cs_end):
    reads_seen_window += 1

    if read.is_unmapped:
        continue
    if read.mapping_quality < MIN_MAPQ:
        continue
    reads_pass_mapq += 1

    if not full_span(read, cs_start, cs_end):
        continue
    reads_pass_fullspan += 1

    seq = read.query_sequence
    pairs = read.get_aligned_pairs(matches_only=False)

    region = []
    has_deletion = False

    for qpos, rpos in pairs:
        if rpos is not None and (cs_start - 1) <= rpos <= (cs_end - 1):
            if qpos is None:
                region.append("-")
                has_deletion = True
            else:
                region.append(seq[qpos])

    if has_deletion:
        reads_with_deletion += 1

    nt = "".join(region)
    aa = translate(nt)
    motifs.append(aa)

bam.close()

# -----------------------------
# Aggregate
# -----------------------------
if len(motifs) == 0:
    counts_all = pd.DataFrame(columns=["Motif", "Reads", "Freq"])
    counts = counts_all.copy()
    total = 0
else:
    df = pd.DataFrame({"Motif": motifs})
    counts_all = (
        df.value_counts()
          .reset_index(name="Reads")
          .sort_values("Reads", ascending=False)
          .reset_index(drop=True)
    )
    total = int(counts_all["Reads"].sum())
    counts_all["Freq"] = counts_all["Reads"] / total
    counts = counts_all[counts_all["Freq"] >= MIN_FREQ].copy()

counts.to_csv(f"{sample}.cleavage_haplotypes.tsv", sep="\\t", index=False)
counts_all.to_csv(f"{sample}.cleavage_haplotypes.all.tsv", sep="\\t", index=False)

summary = pd.DataFrame([
    ["sample", sample],
    ["contig", contig],
    ["cs_start", cs_start],
    ["cs_end", cs_end],
    ["min_mapq", MIN_MAPQ],
    ["min_freq", MIN_FREQ],
    ["reads_seen_window", reads_seen_window],
    ["reads_pass_mapq", reads_pass_mapq],
    ["reads_pass_fullspan", reads_pass_fullspan],
    ["reads_with_deletion_in_cs", reads_with_deletion],
    ["total_motif_reads", total],
    ["unique_motifs_all", len(counts_all)],
    ["reported_motifs_ge_threshold", len(counts)]
], columns=["metric", "value"])

summary.to_csv(f"{sample}.cleavage_haplotypes.summary.tsv", sep="\\t", index=False)
EOF
    """
}