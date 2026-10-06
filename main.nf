#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --- params ---
params.sample = ''
params.outdir = ''

params.fastq_dir = ''
params.reference = ''
params.script_rmdup = "${projectDir}/samremovedup.py"
params.chromosomes   = "${projectDir}/chromosomes.tsv"

// --- variant caller ---
// caller = 'bcftools' : bcftools mpileup/call, sites selected by read counts (default)
// caller = 'dv'       : DeepVariant — runs in container
params.caller = 'bcftools'

// --- two-phase strategy ---
// Phase 1 (SNV-only): converge the reference on substitutions with stable coordinates.
//                     Improves read mapping before tackling indels.
// Phase 2 (SNPs + indels): fix insertions/deletions once the SNP landscape is resolved
//                          and the caller has a cleaner reference to work with.
params.phase1_iterations = 3    // iterations 1–3 = SNV-only; remaining = SNPs + indels

// --- mapping ---
// bwa mem options. Keep the default gap-open (-O 6) and clipping (-L 5) penalties:
// -O 4 -L 8 put 1-bp gaps at read ends (12% of reads at iteration 1), which hides
// the alt allele and leaves differences unfixed.
params.bwa_mem_args = '-k 19 -r 2.5'

// --- variant calling parameters ---
params.mapQ = 25
params.baseQ = 30
params.min_readlength = 40
params.min_depth = 3
params.max_depth = 50
params.variant_quality = 0

// --- bcftools site selection (caller = 'bcftools') ---
// a site is applied to the consensus when it has >= min_alt_reads alt reads
// and alt reads / depth >= min_alt_frac (1.0 = every read must show the alt allele)
params.min_alt_reads = 3
params.min_alt_frac = 0.8

// GATK3 indel realignment before bcftools calling, run in phase 2 (indel
// iterations) only. Without it ~20% fewer indels pass min_alt_frac.
params.indel_realign = true

// bcftools mpileup indel model in phase 2: 'default', 'indels-2.0' or 'indels-cns'
params.indel_model = 'indels-cns'


include { iterative_mapping as iterative_mapping_1 } from './iterative_mapping_workflow.nf'
include { iterative_mapping as iterative_mapping_2 } from './iterative_mapping_workflow.nf'
include { iterative_mapping as iterative_mapping_3 } from './iterative_mapping_workflow.nf'
include { iterative_mapping as iterative_mapping_4 } from './iterative_mapping_workflow.nf'
include { iterative_mapping as iterative_mapping_5 } from './iterative_mapping_workflow.nf'
include { iterative_mapping as iterative_mapping_6 } from './iterative_mapping_workflow.nf'
include { iterative_mapping as iterative_mapping_7 } from './iterative_mapping_workflow.nf'
include { iterative_mapping as iterative_mapping_8 } from './iterative_mapping_workflow.nf'

workflow {
    // =====================================================================
    // Phase 1 : SNV-only iterations (converge substitutions, stable coords)
    // =====================================================================
    ch_iter1 = channel.value(1)
    ch_ref_1 = channel.fromPath(params.reference)
    fastq_files_1 = channel.fromPath("${params.fastq_dir}*.{fastq,fq}.gz")
    iterative_mapping_1(ch_ref_1, fastq_files_1, ch_iter1, channel.value(true))

    ch_iter2 = channel.value(2)
    ch_ref_2 = iterative_mapping_1.out
    fastq_files_2 = channel.fromPath("${params.fastq_dir}*.{fastq,fq}.gz")
    iterative_mapping_2(ch_ref_2, fastq_files_2, ch_iter2, channel.value(true))

    ch_iter3 = channel.value(3)
    ch_ref_3 = iterative_mapping_2.out
    fastq_files_3 = channel.fromPath("${params.fastq_dir}*.{fastq,fq}.gz")
    iterative_mapping_3(ch_ref_3, fastq_files_3, ch_iter3, channel.value(true))

    // =====================================================================
    // Phase 2 : SNPs + indels (fix indels now that SNPs are converged)
    // =====================================================================
    ch_iter4 = channel.value(4)
    ch_ref_4 = iterative_mapping_3.out
    fastq_files_4 = channel.fromPath("${params.fastq_dir}*.{fastq,fq}.gz")
    iterative_mapping_4(ch_ref_4, fastq_files_4, ch_iter4, channel.value(false))

    ch_iter5 = channel.value(5)
    ch_ref_5 = iterative_mapping_4.out
    fastq_files_5 = channel.fromPath("${params.fastq_dir}*.{fastq,fq}.gz")
    iterative_mapping_5(ch_ref_5, fastq_files_5, ch_iter5, channel.value(false))

    // sixth iteration (uncomment to add more indel iterations)
    ch_iter6 = channel.value(6)
    ch_ref_6 = iterative_mapping_5.out
    fastq_files_6 = channel.fromPath("${params.fastq_dir}*.{fastq,fq}.gz")
    iterative_mapping_6(ch_ref_6, fastq_files_6, ch_iter6, channel.value(false))

    // seventh iteration (uncomment to add more indel iterations)
    ch_iter7 = channel.value(7)
    ch_ref_7 = iterative_mapping_6.out
    fastq_files_7 = channel.fromPath("${params.fastq_dir}*.{fastq,fq}.gz")
    iterative_mapping_7(ch_ref_7, fastq_files_7, ch_iter7, channel.value(false))

    ch_iter8 = channel.value(8)
    ch_ref_8 = iterative_mapping_7.out
    fastq_files_8 = channel.fromPath("${params.fastq_dir}*.{fastq,fq}.gz")
    iterative_mapping_8(ch_ref_8, fastq_files_8, ch_iter8, channel.value(false))
}


