# Cohort QC Reporting Workflow (`MAKE_REPORTS`)

The `MAKE_REPORTS` workflow (`workflows/make_reports.nf`) runs when `--make_reports true`. It aggregates per-sample outputs across a single project folder (`--project_qc_folder`) or multiple batch folders (`--alignment_qc_inputs` / `--variant_qc_inputs`), runs cohort-wide checks, evaluates QC flags, and generates a unified interactive HTML report.

## Pipeline Steps

1. **Manifest aggregation**
   - Collects and merges `published_files_manifest/*.csv` files across all input alignment and variant QC directories into unified per-modality TSV manifests (`*_inputs.tsv`), which are also published under `qc_reports/input_files/`.
2. **Pedigree preparation (`MAKE_PED`)**
   - Uses `--ped_file` if provided; otherwise builds `<project_id>.ped` automatically from the collected Somalier manifest, preserving declared sex while treating samples as unrelated.
3. **Cohort relatedness, sex & ancestry (`SOMALIER_RELATE`, `SOMALIER_ANCESTRY`)**
   - `SOMALIER_RELATE` evaluates pairwise IBD/relatedness, duplicate sample detection, and genetic vs. declared sex concordance across all collected `.somalier` profiles.
   - `SOMALIER_ANCESTRY` projects cohort samples onto the 1000 Genomes reference panel (`--somalier_1kgdata`, `--somalier_labels`) to infer continental ancestry.
4. **Alignment QC summary (`SUMMARISE_ALIGN_QC`)**
   - Runs `bin/summarise_align_qc.py` across `mosdepth`, `samtools stats`/`idxstats`, `VerifyBamID2`, `qrs`, and `Somalier` outputs using `--report_target_depths`, `--report_target_baseq`, and `--alignment_qc_criteria`.
   - Emits per-sample JSON, 2-header TSV, binary `qc_flags.tsv`, `qc_criteria.tsv`, and a `datavzrd` config snippet.
5. **Variant QC summary (`SUMMARISE_VARIANT_QC`)**
   - Runs `bin/summarise_var_qc.py` across `bcftools stats` (global and regions), `bcftools idxstats`, and `CHARR` outputs using `--report_var_qual`, `--report_var_dp`, and `--variant_qc_criteria`.
   - Emits per-sample JSON, 2-header TSV, binary `qc_flags.tsv`, `qc_criteria.tsv`, and a `datavzrd` config snippet.
6. **Unified report preparation (`PREPARE_QC_REPORT`)**
   - Runs `bin/create_qc_report.py` to:
     - Evaluate alignment and variant flags jointly and produce a project-level QC summary (`alignment_qc`, `variants_qc`, `total` rows with `N_samples`, `N_PASS`, `N_WARN`, `N_FAIL`, and `failure_rate`).
     - Merge alignment and variant flags into a single `qc_flags.tsv` (reporting `N_warning_flags`, `N_fail_flags`, overall `QC_status`, and all flag columns where at least one sample fails/warns).
     - Generate flag count summaries for stacked Vega-Lite bar plots (separated by `FAIL` and `WARN` severities).
     - Assemble the unified `datavzrd` YAML configuration and custom HTML renderer (`project_qc_summary.js`).
7. **Interactive HTML report (`DATAVZRD`)**
   - Renders the static, self-contained HTML report under `qc_reports/report/`.
   - Views are ordered cleanly in the report dropdown and linked by sample ID:
     1. **Project QC Summary** (default view, styled summary table)
     2. **QC Flags Summary** (separate `FAIL` and `WARN` Vega-Lite bar plots)
     3. **QC Flags** (merged per-sample flag matrix with categorical `PASS`/`WARN`/`FAIL` badges)
     4. **Align QC Metrics**
     5. **Var QC Metrics**
     6. **Align QC Criteria**
     7. **Var QC Criteria**
