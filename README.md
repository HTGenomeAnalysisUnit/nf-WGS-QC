# nf-WGS-QC

**nf-WGS-QC** is a Nextflow pipeline for comprehensive, scalable Whole Genome Sequencing (WGS) quality control across alignments (BAM/CRAM), variant calls (VCF/gVCF/BCF), and cohort-level interactive reporting.

## Features

- **Alignment QC (`ALIGNMENTS_QC`)**
  - Automatic BAM/CRAM indexing (`samtools index`) and sample manifest extraction (`BAM_MANIFEST`).
  - Whole-genome, per-region, threshold, and quantized coverage profiling (`mosdepth`).
  - Alignment summary and per-chromosome read metrics (`samtools stats`, `samtools idxstats`).
  - Per-flowcell/lane read statistics (`Quick Reads Stats` / `qrs`).
  - Sample-level contamination estimation (`VerifyBamID2`).
  - Genetic fingerprint extraction for downstream relatedness, sex, and ancestry checks (`Somalier extract`).
- **Variant QC (`VARIANTS_QC`)**
  - Automatic VCF/gVCF/BCF indexing (`bcftools index`), variant filtering (`bcftools view`), and manifest extraction (`VCF_MANIFEST`).
  - Whole-genome and optional region-restricted variant statistics (`bcftools stats`, `bcftools idxstats`).
  - Variant-level contamination estimation with CHARR (`sceVCF`), including optional allele-frequency annotation via `Echtvar`.
- **Cohort QC Reporting (`MAKE_REPORTS`)**
  - Aggregates outputs across a single project folder or multiple incremental QC batches.
  - Cohort-wide relatedness, duplicate detection, sex concordance (`Somalier relate`), and 1000 Genomes ancestry projection (`Somalier ancestry`).
  - Customizable YAML `FAIL`/`WARN` criteria evaluated across both alignment and variant metrics.
  - Unified, interactive HTML report built with `datavzrd` featuring project QC summaries, flag breakdowns, and sample-linked metric tables.
- **Reproducible & Validated**
  - Parameter and samplesheet validation via `nf-schema`.
  - Containerized execution with Singularity, Apptainer, or Docker.

## Requirements

- [Nextflow](https://www.nextflow.io/) `>= 26.04`
- Container runtime: [Singularity](https://sylabs.io/singularity/), [Apptainer](https://apptainer.org/), or [Docker](https://www.docker.com/)

## Quick Start

### 1. Alignment QC

```bash
nextflow run HTGenomeAnalysisUnit/nf-WGS-QC -profile singularity \
    --project_id my_project \
    --outdir QC_results \
    --genome_fasta /path/to/genome.fasta \
    --input_bam samples_bam.tsv \
    --somalier_vars /path/to/sites.hg38.vcf.gz \
    --verifybamid_svd_file_prefix /path/to/1000g.phase3.100k.b38.vcf.gz.dat
```

### 2. Variant QC

```bash
nextflow run HTGenomeAnalysisUnit/nf-WGS-QC -profile singularity \
    --project_id my_project \
    --outdir QC_results \
    --genome_fasta /path/to/genome.fasta \
    --input_vcf samples_vcf.tsv \
    --af_annotation_zip /path/to/gnomad.echtvar.zip
```

> Alignment QC and Variant QC can be run together in a single command by providing both `--input_bam` and `--input_vcf`.

### 3. Cohort QC Reports

Run report generation separately after Alignment QC and/or Variant QC have completed:

```bash
nextflow run HTGenomeAnalysisUnit/nf-WGS-QC -profile singularity \
    --project_id my_project \
    --outdir QC_results \
    --make_reports true \
    --project_qc_folder QC_results \
    --alignment_qc_criteria assets/align_qc_flags.yml \
    --variant_qc_criteria assets/var_qc_flags.yml \
    --somalier_labels /path/to/ancestry-labels-1kg.tsv \
    --somalier_1kgdata /path/to/1kg-somalier/
```

For full usage details, samplesheet formats, parameter references, and QC rule customization, see the **[nf-WGS-QC Documentation](https://htgenomeanalysisunit.github.io/nf-WGS-QC/)** (source in [`docs/`](docs/)).
