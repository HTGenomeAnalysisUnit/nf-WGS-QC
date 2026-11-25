# Inputs & Samplesheets

All pipeline inputs are validated at startup using `nf-schema` JSON schemas located in `assets/`.

## Alignment QC Samplesheet (`--input_bam`)

Provide a TSV or CSV file describing input BAM or CRAM files.

| Column | Required | Description |
| :--- | :---: | :--- |
| `sample` | Yes | Unique sample identifier (no whitespace). |
| `alignment_file` | Yes | Path to an existing `.bam` or `.cram` file (unique per row). Indexes (`.bai`/`.crai`) are detected automatically if present alongside the file; otherwise `SAMTOOLS_INDEX` creates them. |
| `sex` | No | Declared sample sex: `M`, `F`, `NA`, `1` (male), `2` (female), or `0` (unknown). Default: `0`. |
| `batch_label` | No | Batch identifier string (no whitespace). Default: `project_batch`. |
| `regions_bed` | No | Per-sample regions `.bed` or `.bed.gz` file. Overrides `--regions_bed` when provided. |

Example `samples_bam.tsv`:

```tsv
sample	sex	batch_label	alignment_file	regions_bed
SAMPLE_01	M	batch1	/data/bams/SAMPLE_01.cram	
SAMPLE_02	F	batch1	/data/bams/SAMPLE_02.cram	/data/beds/capture.bed.gz
```

---

## Variant QC Samplesheet (`--input_vcf`)

Provide a TSV or CSV file describing input VCF, gVCF, or BCF files.

| Column | Required | Description |
| :--- | :---: | :--- |
| `sample` | Yes | Unique sample identifier (no whitespace). |
| `variant_file` | Yes | Path to an existing `.vcf`, `.vcf.gz`, `.g.vcf`, `.g.vcf.gz`, `.bcf`, or `.bcf.gz` file (unique per row). Indexes (`.tbi`/`.csi`) are detected automatically or created via `BCFTOOLS_INDEX`. |
| `sex` | No | Declared sample sex (`M`, `F`, `NA`, `1`, `2`, `0`). Default: `0`. |
| `batch_label` | No | Batch identifier string (no whitespace). Default: `project_batch`. |
| `regions_bed` | No | Per-sample regions `.bed` or `.bed.gz` file. Overrides `--regions_bed` when provided. |

Example `samples_vcf.tsv`:

```tsv
sample	sex	batch_label	variant_file
SAMPLE_01	M	batch1	/data/vcfs/SAMPLE_01.g.vcf.gz
SAMPLE_02	F	batch1	/data/vcfs/SAMPLE_02.g.vcf.gz
```

---

## Report Generation Inputs (`--make_reports true`)

Cohort report generation collects published manifests produced by previous `ALIGNMENTS_QC` and `VARIANTS_QC` runs. You can supply inputs in two ways:

1. **Single project folder (`--project_qc_folder`)**: Points to a directory containing `qc_alignments/` and `qc_variants/` subfolders (the default layout produced by `--outdir`).
2. **Multi-batch folder lists (`--alignment_qc_inputs`, `--variant_qc_inputs`)**: Text/TSV/CSV files with one existing `qc_alignments` or `qc_variants` directory path per line (no header). Useful for combining batches processed separately over time.

Example `alignment_qc_folders.txt`:

```text
/projects/cohortA/batch1/qc_alignments
/projects/cohortA/batch2/qc_alignments
```

---

## Reference & Auxiliary Files

| Parameter | Workflow(s) | Requirement |
| :--- | :--- | :--- |
| `--genome_fasta` | Alignment QC, Variant QC | Reference genome `.fasta` with a matching `.fasta.fai` index in the same directory. |
| `--regions_bed` | Alignment QC, Variant QC | Optional `.bed` / `.bed.gz` used for region-restricted coverage and variant stats when not specified per sample in the samplesheet. |
| `--somalier_vars` | Alignment QC | Somalier polymorphic sites VCF (e.g., `sites.hg38.vcf.gz`). |
| `--verifybamid_svd_file_prefix` | Alignment QC | Prefix ending in `.dat` for VerifyBamID2 SVD resource files (`.UD`, `.mu`, `.bed` must exist at `<prefix>.UD`, `<prefix>.mu`, `<prefix>.bed`). |
| `--af_annotation_zip` | Variant QC | Echtvar `.zip` allele-frequency archive (required when `--annotate_af true`, the default). |
| `--somalier_labels` | Make Reports | TSV file mapping 1000 Genomes reference sample IDs to ancestry labels (`ancestry-labels-1kg.tsv`). |
| `--somalier_1kgdata` | Make Reports | Directory containing `.somalier` reference extracts for 1000 Genomes ancestry projection. |
| `--ped_file` | Make Reports | Optional PED file describing expected family relationships. If omitted, `MAKE_PED` generates a PED file treating all samples as unrelated while preserving declared sex. |
| `--alignment_qc_criteria` | Make Reports | YAML file defining alignment QC `fail_rules` and `warning_rules` (see [assets/align_qc_flags.yml](https://github.com/HTGenomeAnalysisUnit/nf-WGS-QC/blob/main/assets/align_qc_flags.yml)). |
| `--variant_qc_criteria` | Make Reports | YAML file defining variant QC `fail_rules` and `warning_rules` (see [assets/var_qc_flags.yml](https://github.com/HTGenomeAnalysisUnit/nf-WGS-QC/blob/main/assets/var_qc_flags.yml)). |
