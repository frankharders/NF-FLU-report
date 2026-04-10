WBVR specific report from NF-Flu pipeline

no version yet

all used modules are in a separate folder
All started from main.

command:
nextflow run main.nf -resume   --samplesheet  /mnt/lely_archive/harde004/backup/backup-2026/AI-2026/NGS-26-AP/samplesheet.csv   --genin2       /mnt/lely_archive/harde004/backup/backup-2026/AI-2026/NGS-26-AP/results/genin2/genin2.tsv   --annotdir     /mnt/lely_archive/harde004/backup/backup-2026/AI-2026/NGS-26-AP/results/annotation   --flumut       /mnt/lely_archive/harde004/backup/backup-2026/AI-2026/NGS-26-AP/results/flumut/flumut-markers.tsv   --consensusdir /mnt/lely_archive/harde004/backup/backup-2026/AI-2026/NGS-26-AP/results/consensus/bcftools   --outdir       /mnt/lely_archive/harde004/backup/backup-2026/AI-2026/NGS-26-AP/results/report

results for the report are in 

tree results/report/
results/report/
├── 2026-wk03-01_mutation_check.tsv
├── 2026-wk03-02_mutation_check.tsv
├── bam
│   ├── 2026-wk03-01.sorted.bam
│   ├── 2026-wk03-01.sorted.bam.bai
│   ├── 2026-wk03-02.sorted.bam
│   └── 2026-wk03-02.sorted.bam.bai
├── cleavage_site
│   ├── 2026-wk03-01.cleavage_haplotypes.all.tsv
│   ├── 2026-wk03-01.cleavage_haplotypes.summary.tsv
│   ├── 2026-wk03-01.cleavage_haplotypes.tsv
│   ├── 2026-wk03-02.cleavage_haplotypes.all.tsv
│   ├── 2026-wk03-02.cleavage_haplotypes.summary.tsv
│   └── 2026-wk03-02.cleavage_haplotypes.tsv
├── consensus_renamed
│   ├── 2026-wk03-01.renamed.consensus.fasta
│   └── 2026-wk03-02.renamed.consensus.fasta
├── final_report.tsv
├── flumut_parsed.tsv
├── plots
│   ├── 2026-wk03-01_coverage.pdf
│   └── 2026-wk03-02_coverage.pdf
├── simple_report.tsv
├── trimmed
│   ├── 2026-wk03-01_fastp.html
│   ├── 2026-wk03-01_fastp.json
│   ├── 2026-wk03-01_trimmed_R1.fastq.gz
│   ├── 2026-wk03-01_trimmed_R2.fastq.gz
│   ├── 2026-wk03-02_fastp.html
│   ├── 2026-wk03-02_fastp.json
│   ├── 2026-wk03-02_trimmed_R1.fastq.gz
│   ├── 2026-wk03-02_trimmed_R2.fastq.gz
│   └── test
│       ├── 2026-wk03-01_trimmed_R1.fastq.gz
│       └── R1.fq.gz
└── variants
    ├── 2026-wk03-01.depth.tsv
    ├── 2026-wk03-01.vcf.gz
    ├── 2026-wk03-01.vcf.gz.tbi
    ├── 2026-wk03-02.depth.tsv
    ├── 2026-wk03-02.vcf.gz
    └── 2026-wk03-02.vcf.gz.tbi


collaboration:
Marc Engelsma, Rene Heutink & Frank Harders

# end 10 April 2026

