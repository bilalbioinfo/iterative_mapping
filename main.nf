#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --- params ---
params.sample = ''
params.outdir = ''

params.fastq_dir = ''
params.reference = ''
params.script_rmdup = './samremovedup.py'
params.chromosomes   = './chromosomes.tsv'

// --- variant caller ---
// caller = 'gatk' : GATK HaplotypeCaller (default)
// caller = 'dv'   : DeepVariant — faster, runs in container
params.caller = 'dv'

// --- two-phase strategy ---
// Phase 1 (SNV-only): converge the reference on substitutions with stable coordinates.
//                     Improves read mapping before tackling indels.
// Phase 2 (SNPs + indels): fix insertions/deletions once the SNP landscape is resolved
//                          and HaplotypeCaller has a cleaner reference to work with.
params.phase1_iterations = 3    // iterations 1–3 = SNV-only; remaining = SNPs + indels

// --- variant calling parameters ---
params.mapQ = 20
params.baseQ = 30
params.min_readlength = 40
params.min_depth = 3
params.max_depth = 50
params.variant_quality = 0


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
    ch_iter1 = Channel.value(1)
    ch_ref_1 = Channel.fromPath(params.reference)
    fastq_files_1 = Channel.fromPath("${params.fastq_dir}*.fastq.gz")
    iterative_mapping_1(ch_ref_1, fastq_files_1, ch_iter1, Channel.value(true))

    ch_iter2 = Channel.value(2)
    ch_ref_2 = iterative_mapping_1.out.consensus_fasta
    fastq_files_2 = Channel.fromPath("${params.fastq_dir}*.fastq.gz")
    iterative_mapping_2(ch_ref_2, fastq_files_2, ch_iter2, Channel.value(true))

    ch_iter3 = Channel.value(3)
    ch_ref_3 = iterative_mapping_2.out.consensus_fasta
    fastq_files_3 = Channel.fromPath("${params.fastq_dir}*.fastq.gz")
    iterative_mapping_3(ch_ref_3, fastq_files_3, ch_iter3, Channel.value(true))

    // =====================================================================
    // Phase 2 : SNPs + indels (fix indels now that SNPs are converged)
    // =====================================================================
    ch_iter4 = Channel.value(4)
    ch_ref_4 = iterative_mapping_3.out.consensus_fasta
    fastq_files_4 = Channel.fromPath("${params.fastq_dir}*.fastq.gz")
    iterative_mapping_4(ch_ref_4, fastq_files_4, ch_iter4, Channel.value(false))

    ch_iter5 = Channel.value(5)
    ch_ref_5 = iterative_mapping_4.out.consensus_fasta
    fastq_files_5 = Channel.fromPath("${params.fastq_dir}*.fastq.gz")
    iterative_mapping_5(ch_ref_5, fastq_files_5, ch_iter5, Channel.value(false))

    // sixth iteration (uncomment to add more indel iterations)
    ch_iter6 = Channel.value(6)
    ch_ref_6 = iterative_mapping_5.out.consensus_fasta
    fastq_files_6 = Channel.fromPath("${params.fastq_dir}*fastq.gz")
    iterative_mapping_6(ch_ref_6, fastq_files_6, ch_iter6, Channel.value(false))

    // seventh iteration (uncomment to add more indel iterations)
    ch_iter7 = Channel.value(7)
    ch_ref_7 = iterative_mapping_6.out.consensus_fasta
    fastq_files_7 = Channel.fromPath("${params.fastq_dir}*fastq.gz")
    iterative_mapping_7(ch_ref_7, fastq_files_7, ch_iter7, Channel.value(false))

    ch_iter8 = Channel.value(8)
    ch_ref_8 = iterative_mapping_7.out.consensus_fasta
    fastq_files_8 = Channel.fromPath("${params.fastq_dir}*fastq.gz")
    iterative_mapping_8(ch_ref_8, fastq_files_8, ch_iter8, Channel.value(false))
}


