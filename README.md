# Iterative Mapping Pipeline

This pipeline is for building a sample-specific consensus reference from short-read sequencing data. The idea is simple: instead of mapping reads once to a distant reference genome, you map, call variants, update the reference with those variants, and repeat. Each round the reference gets closer to the actual sample, so reads map better and variants are called more accurately.

This is useful for ancient DNA or any situation where the sample is diverged enough from the reference that a single round of mapping leaves a lot on the table.

## What it does

The pipeline runs 8 mapping iterations in two phases:

**Phase 1 (iterations 1–3)** — SNVs only. It calls only substitutions and uses those to update the reference. Indels are excluded here on purpose because they shift coordinates and can cause instability early on. Getting the SNP landscape right first makes everything downstream more reliable.

**Phase 2 (iterations 4–8)** — SNPs + indels. Once the reference has converged on substitutions, indels are introduced. By this point the reference is close enough that indel calls are actually trustworthy.

Each iteration does the same thing: index the reference → map all FASTQs with BWA-MEM → filter by mapping quality and read length → remove duplicates → call variants per chromosome → merge variant calls → generate a new consensus FASTA. That consensus becomes the reference for the next iteration.

## Variant callers

We recommend bcftools (the default, `--caller bcftools`). Reads are piled up with `bcftools mpileup -B` and a site is applied to the consensus when it has at least `--min_alt_reads` reads showing the alt allele and those reads make up at least `--min_alt_frac` of the depth. The decision is based on read counts only, not on a genotype model.

BAQ is switched off on purpose. In a sample that is strongly diverged from the reference, mismatches sit close together, and BAQ lowers the quality of exactly those bases. In a test on Neanderthal reads mapped to macaque, BAQ dropped about 5% of the sites where every read showed the alt allele.

In phase 2, reads are realigned around indels with GATK3 IndelRealigner before calling (`--indel_realign`, on by default). Without realignment, reads whose indel lies near the read end are aligned with mismatches instead of a gap, so the indel never reaches `--min_alt_frac`. In the test this cost about 20% of the indels. SNV calls are almost unchanged by realignment (under 1% of sites differ), so it is skipped in phase 1. Realignment takes 30–60 minutes per chromosome per iteration, so you can turn it off with `--indel_realign false` to save time at the cost of fewer indels.

DeepVariant (`--caller dv`) is also supported. Its model is trained on human data and is cautious with reads that carry many close mismatches. In the same test it returned many sites as `./.` (NoCall) even when every read showed the alt allele, so those sites were never applied. It runs in an Apptainer container.

## Requirements

- Nextflow 26.x (strict syntax). Also runs on 25.x
- SLURM cluster (configured for NAISS/Dardel by default — change the partition and account in `nextflow.config`)
- BWA, SAMtools, bcftools available as modules
- Apptainer with the GATK3 image (`broadinstitute/gatk3:3.8-1`) for indel realignment, and the DeepVariant image for `--caller dv`
- Python 3 (for the duplicate removal script `samremovedup.py`)

## Input files

- FASTQ files: one or more `*.fastq.gz` files in a single directory
- Reference genome: a FASTA file to start from
- Chromosomes list: `chromosomes.tsv` — one chromosome/scaffold name per line, must match the sequence names in your reference FASTA. Lines starting with `#` are ignored

The chromosomes file controls which regions get variant calls. Edit it to match your reference. There is an example in the file for domestic cat.

## How to run

```bash
nextflow run main.nf \
  --sample SAMPLE_NAME \
  --outdir /path/to/output/ \
  --fastq_dir /path/to/fastqs/ \
  --reference /path/to/reference.fasta \
  --chromosomes chromosomes.tsv \
  --script_rmdup /path/to/samremovedup.py
```

Resume a run after failure or interruption:

```bash
nextflow run main.nf -resume [same parameters]
```

Resume is enabled by default in `nextflow.config`.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--sample` | (required) | Sample name, used in output file names |
| `--outdir` | (required) | Output directory |
| `--fastq_dir` | (required) | Directory containing input FASTQ files (`*.fastq.gz`) |
| `--reference` | (required) | Starting reference FASTA |
| `--caller` | `bcftools` | Variant caller: `bcftools` or `dv` (DeepVariant) |
| `--phase1_iterations` | `3` | How many SNV-only iterations before switching to SNP+indel mode |
| `--bwa_mem_args` | `-k 19 -r 2.5` | Options passed to `bwa mem`. Avoid lowering the gap-open (`-O`) or raising the clipping (`-L`) penalty: that creates 1-bp gaps at read ends |
| `--mapQ` | `25` | Minimum mapping quality |
| `--baseQ` | `30` | Minimum base quality (bcftools only) |
| `--min_readlength` | `40` | Minimum read length to keep after mapping |
| `--min_alt_reads` | `3` | Minimum reads showing the alt allele (bcftools only) |
| `--min_alt_frac` | `0.8` | Minimum fraction of reads showing the alt allele (with 3 reads all must agree); `1.0` = every read (bcftools only) |
| `--indel_realign` | `true` | GATK3 indel realignment before bcftools calling in phase 2 (bcftools only) |
| `--indel_model` | `indels-cns` | bcftools mpileup indel model in phase 2: `default`, `indels-2.0` or `indels-cns` (bcftools only) |
| `--min_depth` | `3` | Minimum depth to call a variant (dv only) |
| `--max_depth` | `50` | Maximum depth to call a variant (caps repetitive regions) |

## Output structure

Each iteration produces its own directory under `--outdir`:

```
outdir/
  iteration_1/
    raw_bams/        per-library sorted, filtered, deduplicated BAMs (symlinks into work/)
    merged_bams/     all libraries merged into one BAM (symlinks into work/)
    variants/        per-chromosome BCFs and the merged genome-wide BCF (copied)
    reference/       the consensus FASTA used as input for the next iteration (copied)
  iteration_2/
    ...
  iteration_8/
    reference/       final consensus — this is your end result
```

The final consensus FASTA from iteration 8 is your sample-specific reference. Consensus FASTAs, logs and variant BCFs are copied into `--outdir`, so `work/` can be deleted once the run has finished; the BAM links then break, but BAMs are intermediates (remap reads to the final consensus for downstream work). You can use that for any downstream analysis that benefits from a closer reference.

## Cluster configuration

The pipeline is configured to run on SLURM. The partition and account are set in `nextflow.config` — change them to match your allocation. Resource limits per process (CPUs, memory, time) are also defined there and can be adjusted if jobs are timing out or getting killed.

Limits are set for a whole mammalian genome (~2.5 Gb, ~35x pooled coverage): per-chromosome steps get up to 12 h, merging and duplicate removal of all reads up to 24 h. DeepVariant jobs are the most resource-heavy (32 CPUs, 64 GB RAM) because they run per chromosome. If your reference has many chromosomes, consider batching or increasing the queue size.
