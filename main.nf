#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PARSE_INPUTS       } from './modules/parse_inputs'
include { MAKE_REPORT        } from './modules/make_report'
include { PARSE_FLUMUT       } from './modules/parse_flumut'
include { CHECK_MUTATIONS    } from './modules/check_mutations'
include { DETECT_PLATFORM    } from './modules/detect_platform'
include { TRIM_READS         } from './modules/trim_reads'
include { MAP_READS          } from './modules/map_reads'
include { CALL_VARIANTS      } from './modules/call_variants'
include { COVERAGE_PLOT      } from './modules/coverage_plot'
include { CHECK_CLEAVAGESITE } from './modules/check_cleavagesite'
include { EXPORT_CONSENSUS   } from './modules/export_consensus'

params.samplesheet  = params.samplesheet  ?: "$PWD/samplesheet.csv"
params.genin2       = params.genin2       ?: "$PWD/results/genin2/genin2.tsv"
params.annotdir     = params.annotdir     ?: "$PWD/results/annotation"
params.flumut       = params.flumut       ?: "$PWD/results/flumut/flumut-markers.tsv"
params.consensusdir = params.consensusdir ?: "$PWD/results/consensus/bcftools"
params.outdir       = params.outdir       ?: "$PWD/results/report"

params.cs_min_mapq  = params.cs_min_mapq  ?: 20
params.cs_min_freq  = params.cs_min_freq  ?: 0.01

workflow {

    ch_samplesheet = Channel.fromPath(params.samplesheet)

    /*
     * Platform detectie
     */
    DETECT_PLATFORM(ch_samplesheet)
    DETECT_PLATFORM.out.platform.view { p -> "Platform: $p" }

    /*
     * Reads kanaal
     */
    ch_reads = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (row.containsKey('fastq_2') && row.fastq_2) {
                tuple(row.sample, 'illumina', [file(row.fastq_1), file(row.fastq_2)])
            } else {
                tuple(row.sample, 'nanopore', [file(row.reads)])
            }
        }

    TRIM_READS(ch_reads)

    /*
     * Consensus per sample
     */
    ch_consensus = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def consensus = file("${params.consensusdir}/${row.sample}.consensus.fasta")
            tuple(row.sample, consensus)
        }

    /*
     * Export renamed consensus fasta
     */
    EXPORT_CONSENSUS(ch_consensus)

    /*
     * Map reads
     */
    ch_map_input = TRIM_READS.out.trimmed_reads
        .join(ch_consensus, by: 0)
        .map { sampleName, platform, reads, consensus ->
            tuple(tuple(sampleName, platform, reads), consensus)
        }

    MAP_READS(
        ch_map_input.map { it[0] },
        ch_map_input.map { it[1] }
    )

    /*
     * Variant calling
     */
    ch_variants_input = MAP_READS.out.bam
        .join(ch_consensus, by: 0)
        .map { sampleName, bam, bai, consensus ->
            tuple(tuple(sampleName, bam, bai), consensus)
        }

    CALL_VARIANTS(
        ch_variants_input.map { it[0] },
        ch_variants_input.map { it[1] }
    )

    /*
     * Read-level cleavage-site analyse
     * Gebruikt sample-specifieke GFF:
     * ${params.annotdir}/bcftools/<sample>/<sample>.gff
     */
    ch_cleavage_gff = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def gff = file("${params.annotdir}/bcftools/${row.sample}/${row.sample}.gff")
            tuple(row.sample, gff)
        }

    ch_check_cleavage_input = MAP_READS.out.bam
        .join(ch_cleavage_gff, by: 0)
        .map { sampleName, bam, bai, gff ->
            tuple(sampleName, bam, bai, gff)
        }

    CHECK_CLEAVAGESITE(ch_check_cleavage_input)

    /*
     * Bestaande modules
     */
    ch_samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def cleavage = file("${params.annotdir}/${row.sample}/${row.sample}.cleavage.tsv")
            tuple(row.sample, cleavage)
        }

    ch_genin2 = Channel.fromPath(params.genin2)
    ch_flumut = Channel.fromPath(params.flumut)
    ch_input  = ch_samples.combine(ch_genin2)

    PARSE_INPUTS(ch_input)
    PARSE_FLUMUT(ch_flumut)

    ch_faa = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def faa = file("${params.annotdir}/bcftools/${row.sample}/${row.sample}.faa")
            tuple(row.sample, faa)
        }

    ch_check = ch_faa.combine(PARSE_FLUMUT.out.flumut_parsed)
    CHECK_MUTATIONS(ch_check)

    MAKE_REPORT(
        PARSE_INPUTS.out.report_input.collect(),
        CHECK_MUTATIONS.out.mutation_check.collect()
    )

    /*
     * mutation_check met sampleName voor join
     */
    ch_mutation_check = CHECK_MUTATIONS.out.mutation_check
        .map { f ->
            def sampleName = f.name.replace("_mutation_check.tsv", "")
            tuple(sampleName, f)
        }

    /*
     * Plot script als input
     */
    ch_plot_script = Channel.fromPath("$projectDir/modules/coverage_plot.py")

    /*
     * Coverage plot input
     */
    CALL_VARIANTS.out.vcf
        .join(CALL_VARIANTS.out.depth, by: 0)
        .map { sampleName, vcf, tbi, depth ->
            def gff = file("${params.annotdir}/bcftools/${sampleName}/${sampleName}.gff")
            tuple(sampleName, vcf, tbi, depth, gff)
        }
        .join(ch_mutation_check, by: 0)
        .map { sampleName, vcf, tbi, depth, gff, mutation_check ->
            tuple(
                tuple(sampleName, vcf, tbi),
                tuple(sampleName, depth),
                gff,
                mutation_check
            )
        }
        .combine(ch_plot_script)
        .map { vcf_tuple, depth_tuple, gff, mutation_check, plot_script ->
            tuple(vcf_tuple, depth_tuple, gff, mutation_check, plot_script)
        }
        .set { ch_coverage_input }

    COVERAGE_PLOT(
        ch_coverage_input.map { it[0] },
        ch_coverage_input.map { it[1] },
        ch_coverage_input.map { it[2] },
        ch_coverage_input.map { it[3] },
        ch_coverage_input.map { it[4] }
    )
}