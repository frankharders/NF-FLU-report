import re
import csv
import sys
import pysam
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
sample   = sys.argv[1]
vcf_file = sys.argv[2]
dep_file = sys.argv[3]
gff_file = sys.argv[4]
mut_file = sys.argv[5]
depth_df = pd.read_csv(dep_file, sep="\t", header=None,
                        names=["segment", "pos", "depth"])
depth_df = depth_df.dropna()
depth_df["pos"]   = depth_df["pos"].astype(int)
depth_df["depth"] = depth_df["depth"].astype(int)
def seg_order(x):
    m = re.search(r'segment(\d+)', x)
    return int(m.group(1)) if m else 99
segments = sorted(depth_df["segment"].unique().tolist(), key=seg_order)
variants = {}
vf = pysam.VariantFile(vcf_file)
for rec in vf.fetch():
    seg = rec.chrom
    pos = rec.pos
    af  = rec.info.get("AF", [0])[0] if "AF" in rec.info else 0
    variants.setdefault(seg, []).append((pos, float(af)))
cleavage_regions = {}
cds_coords       = {}
with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 9:
            continue
        seg   = cols[0]
        feat  = cols[2]
        start = int(cols[3])
        end   = int(cols[4])
        attrs = cols[8]
        if feat == "misc_feature" and "cleavage site" in attrs:
            cleavage_regions[seg] = (start, end)
        if feat == "CDS":
            m = re.search(r'gene=([^;]+)', attrs)
            if m:
                cds_coords[m.group(1)] = (seg, start)
flumut_positions = {}
with open(mut_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["Sample"] != sample:
            continue
        if row["Mutation_Present"] != "YES":
            continue
        gene = row["Gene"].strip()
        pos  = int(row["Position"].strip())
        if gene not in cds_coords:
            continue
        seg, cds_start = cds_coords[gene]
        nuc_pos = cds_start + (pos - 1) * 3
        flumut_positions.setdefault(seg, []).append(nuc_pos)
antigenic_sites = {
    "Sa" : (163, 185), "Sb"  : (186, 225),
    "Ca1": (169, 204), "Ca2" : (137, 142), "Cb": (70, 75),
}
LOW_COV = 10
n = len(segments)
with PdfPages(f"{sample}_coverage.pdf") as pdf:
    fig, axes = plt.subplots(n, 1, figsize=(18, 4 * n), squeeze=False)
    fig.suptitle(f"Coverage rapport: {sample}", fontsize=16, fontweight="bold")
    for idx, seg in enumerate(segments):
        ax      = axes[idx][0]
        seg_df  = depth_df[depth_df["segment"] == seg].sort_values("pos")
        pos_arr = seg_df["pos"].values.astype(float)
        dep_arr = seg_df["depth"].values.astype(float)
        if len(pos_arr) == 0:
            ax.set_title(seg)
            continue
        mean_cov  = dep_arr.mean()
        seg_label = re.search(r'segment\d+_(\w+)', seg)
        seg_label = seg_label.group(1) if seg_label else seg
        ax.fill_between(pos_arr, dep_arr, alpha=0.4, color="steelblue", label="Coverage")
        ax.plot(pos_arr, dep_arr, color="steelblue", linewidth=0.5)
        ax.axhline(mean_cov, color="navy", linestyle="--", linewidth=1,
                   label=f"Mean: {mean_cov:.0f}x")
        low_mask = dep_arr < LOW_COV
        if low_mask.any():
            ax.fill_between(pos_arr, dep_arr, where=low_mask,
                            alpha=0.4, color="red", label=f"Coverage < {LOW_COV}x")
        snp_labels = set()
        for p, af in variants.get(seg, []):
            if 0.01 <= af < 0.05:
                lbl = "SNP 1-5%" if "SNP 1-5%" not in snp_labels else "_"
                snp_labels.add("SNP 1-5%")
                ax.axvline(p, color="orange", alpha=0.7, linewidth=1.0, label=lbl)
            elif af >= 0.05:
                lbl = "SNP >=5%" if "SNP >=5%" not in snp_labels else "_"
                snp_labels.add("SNP >=5%")
                ax.axvline(p, color="red", alpha=0.9, linewidth=1.5, label=lbl)
        if seg in cleavage_regions:
            cs, ce = cleavage_regions[seg]
            ax.axvspan(cs, ce, alpha=0.25, color="green", label="Cleavage site")
            ax.text((cs+ce)/2, mean_cov * 1.1, "CS",
                    ha="center", fontsize=7, color="darkgreen", fontweight="bold")
        flu_added = False
        for nuc_pos in flumut_positions.get(seg, []):
            lbl = "FluMut YES" if not flu_added else "_"
            flu_added = True
            ax.axvline(nuc_pos, color="blue", alpha=0.8, linewidth=1.5,
                       linestyle=":", label=lbl)
            ax.annotate("▲", xy=(nuc_pos, mean_cov * 1.15),
                        ha="center", fontsize=8, color="blue")
        if seg_label == "HA" and "HA" in cds_coords:
            ha_seg, ha_start = cds_coords["HA"]
            for site, (aa_s, aa_e) in antigenic_sites.items():
                nuc_s = ha_start + (aa_s - 1) * 3
                nuc_e = ha_start + (aa_e - 1) * 3
                ax.axvspan(nuc_s, nuc_e, alpha=0.12, color="purple")
                ax.text((nuc_s+nuc_e)/2, mean_cov * 1.25, site,
                        ha="center", fontsize=6, color="purple", fontweight="bold")
        ax.set_title(f"{seg_label}  |  mean: {mean_cov:.0f}x  |  len: {int(pos_arr[-1])} bp",
                     fontsize=10, fontweight="bold")
        ax.set_xlabel("Positie (nt)", fontsize=8)
        ax.set_ylabel("Coverage (x)", fontsize=8)
        ax.set_xlim(pos_arr[0], pos_arr[-1])
        handles, labels = ax.get_legend_handles_labels()
        seen = {}
        for h, l in zip(handles, labels):
            if l not in seen and not l.startswith("_"):
                seen[l] = h
        ax.legend(seen.values(), seen.keys(), fontsize=7, loc="upper right")
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    pdf.savefig(fig)
    plt.close()
print(f"PDF gegenereerd: {sample}_coverage.pdf")