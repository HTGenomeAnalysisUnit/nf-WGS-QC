# Parameters Reference

All parameters are defined in `nextflow.config` and validated against `nextflow_schema.json` and workflow-specific schemas in `assets/`.

## Output Options

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--outdir` | directory | `QC_results` | Output directory for published results and `pipeline_info/`. |
| `--project_id` | string | *(required)* | Project identifier used in output filenames, manifests, and report headers. |

## Reference Files

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--genome_fasta` | file | *(required for QC)* | Path to reference genome `.fasta` (matching `.fasta.fai` must exist alongside it). |
| `--regions_bed` | file | `null` | Optional `.bed`/`.bed.gz` of regions of interest for region-specific coverage and variant stats. |

## Alignment QC Inputs & Configuration

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--input_bam` | file | `null` | Samplesheet TSV/CSV of BAM/CRAM files (`assets/schema_input_bam.json`). |
| `--somalier_vars` | file | `null` | VCF of Somalier reference polymorphic sites. |
| `--verifybamid_svd_file_prefix` | string | `null` | Prefix (ending in `.dat`) for VerifyBamID2 SVD files (`.UD`, `.mu`, `.bed`). |
| `--mosdepth_precision` | integer | `4` | Decimal precision (`1`–`10`) for `mosdepth` coverage output. |
| `--mosdepth_per_base_coverage` | boolean | `false` | Compute per-base coverage with `mosdepth` (significantly increases runtime). |
| `--mosdepth_thresholds` | string | `"1,10,20,100"` | Comma-separated depth thresholds reported per region when a regions BED is provided. |
| `--mosdepth_quantize_bins` | string | `"0:1:10:20:300:"` | Colon-separated coverage bin boundaries for `mosdepth` quantized output. |
| `--mosdepth_quantize_labels` | string | `"NO_COVERAGE,LOW_COVERAGE,MID_COVERAGE,GOOD_COVERAGE,EXCESSIVE_COVERAGE"` | Labels matching the number of bins in `--mosdepth_quantize_bins`. |
| `--qrs_subset_reads` | integer | `5000000` | Number of reads (`10000`–`10000000`) to subsample in `Quick Reads Stats`. |

## Variant QC Inputs & Configuration

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--input_vcf` | file | `null` | Samplesheet TSV/CSV of VCF/gVCF/BCF files (`assets/schema_input_vcf.json`). |
| `--af_annotation_zip` | file | `null` | Echtvar `.zip` archive used for allele frequency annotation (required when `--annotate_af true`). |
| `--variant_stats_filter` | string | `".,PASS"` | Filter expression passed to `bcftools` (`-f`) when selecting variants for stats. |
| `--annotate_af` | boolean | `true` | Annotate variants with `Echtvar` (after `bcftools norm`) before running CHARR contamination estimation. |
| `--charr_af_field` | string | `"gnomad_af"` | VCF `INFO` field name containing population allele frequencies for `sceVCF`. |

## Report Generation Inputs & Configuration

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--make_reports` | boolean | `false` | Run the `MAKE_REPORTS` cohort aggregation and HTML reporting workflow. |
| `--project_qc_folder` | directory | `null` | Path to a project QC folder containing `qc_alignments/` and `qc_variants/` subfolders. |
| `--alignment_qc_inputs` | file | `null` | Text/TSV/CSV listing one or more `qc_alignments` folders to aggregate. |
| `--variant_qc_inputs` | file | `null` | Text/TSV/CSV listing one or more `qc_variants` folders to aggregate. |
| `--alignment_qc_criteria` | file | `null` | YAML file specifying alignment QC `fail_rules` and `warning_rules`. |
| `--variant_qc_criteria` | file | `null` | YAML file specifying variant QC `fail_rules` and `warning_rules`. |
| `--ped_file` | file | `null` | Optional PED file describing expected sample relationships. |
| `--somalier_labels` | file | `null` | TSV file with 1000 Genomes ancestry labels for `SOMALIER_ANCESTRY`. |
| `--somalier_1kgdata` | directory | `null` | Directory of 1000 Genomes `.somalier` reference extracts. |
| `--report_target_depths` | string | `"10,20,30"` | Comma-separated target depths summarised in the alignment QC report. |
| `--report_target_baseq` | string | `"30"` | Comma-separated base quality thresholds summarised in the alignment QC report. |
| `--report_var_qual` | integer | `20` | Variant `QUAL` threshold summarised in the variant QC report. |
| `--report_var_dp` | integer | `10` | Genotype depth (`DP`) threshold summarised in the variant QC report. |

## Advanced Options

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--help` | boolean | `false` | Print pipeline help message and exit. |
| `--publish_dir_mode` | string | `"copy"` | Nextflow publish mode (`copy`, `move`, `link`, `symlink`, `rellink`). |
| `--trace_report_suffix` | string | `yyyy-MM-dd_HH-mm-ss` | Timestamp suffix for execution timeline, report, trace, DAG, and `pipeline_info/` files. |
| `--monochrome_logs` | boolean | `false` | Disable ANSI colors in pipeline logs. |
| `--custom_config_version` | string | `"master"` | Branch/tag of `nf-core/configs` to load. |
| `--custom_config_base` | string | `https://raw.githubusercontent.com/nf-core/configs/master` | Base URL or local path for institutional configs. |
| `--gpu_container_options` | string | `null` | Extra container runtime options for GPU-enabled processes. |
