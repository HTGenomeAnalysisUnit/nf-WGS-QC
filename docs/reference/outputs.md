# Output Directory Structure

All published results are written under `--outdir` (default: `QC_results/`) using Nextflow's workflow `output` block (`--publish_dir_mode`, default `copy`).

```text
<outdir>/
├── pipeline_info/
│   ├── collated_versions.<timestamp>.yml
│   ├── execution_report_<timestamp>.html
│   ├── execution_timeline_<timestamp>.html
│   ├── execution_trace_<timestamp>.txt
│   ├── params.<timestamp>.yml
│   └── pipeline_dag_<timestamp>.html
│
├── qc_alignments/                          # Published when --input_bam is set
│   ├── input_bam_manifest.tsv
│   ├── sample_to_user_map.tsv
│   ├── coverage/
│   │   ├── summary/                        # mosdepth global/regions dist & summary
│   │   └── detailed/                       # quantized, thresholds, per-base beds & indexes
│   ├── samtools/
│   │   ├── stats/                          # samtools stats outputs
│   │   └── idxstats/                       # samtools idxstats outputs
│   ├── qrs/
│   │   ├── run_files/                      # per_sample_run.tsv, rundetails.tsv
│   │   └── read_stats/                     # read_stats.tsv
│   ├── somalier/
│   │   └── extract/                        # per-sample .somalier profiles
│   ├── verifybamid/                        # .selfSM and .log files
│   └── published_files_manifest/           # Tab-separated manifest indexes for MAKE_REPORTS
│       ├── coverage_summary.csv
│       ├── coverage_detailed.csv
│       ├── samtools_stats.csv
│       ├── samtools_idxstats.csv
│       ├── somalier_extract.csv
│       ├── qrs_run_files.csv
│       ├── qrs_read_stats.csv
│       └── verifybamid_result.csv
│
├── qc_variants/                            # Published when --input_vcf is set
│   ├── input_vcf_manifest.tsv
│   ├── sample_to_user_map.tsv
│   ├── bcftools/
│   │   ├── stats/                          # global & regions bcftools stats
│   │   └── idxstats/                       # per-chromosome variant counts
│   ├── charr_contamination/                # sceVCF CHARR contamination stats
│   └── published_files_manifest/           # Tab-separated manifest indexes for MAKE_REPORTS
│       ├── bcftools_stats.csv
│       ├── bcftools_stats_regions.csv
│       ├── bcftools_idxstats.csv
│       └── charr_contamination.csv
│
└── qc_reports/                             # Published when --make_reports true
    ├── input_files/
    │   ├── <project_id>.ped                # Generated or provided PED file
    │   ├── align_qc/*_inputs.tsv           # Collected alignment QC manifests
    │   └── variant_qc/*_inputs.tsv         # Collected variant QC manifests
    ├── somalier/
    │   ├── relate/                         # pairs.tsv, samples.tsv, HTML report
    │   └── ancestry/                       # somalier-ancestry.tsv, HTML report
    ├── alignments/                         # summarise_align_qc.py TSV, JSON, flags, criteria
    ├── variants/                           # summarise_var_qc.py TSV, JSON, flags, criteria
    └── report/                             # Unified project summary TSVs, datavzrd YAML & HTML report
```

## Key Report Files (`qc_reports/report/`)

- **`<project_id>/index.html`** (inside the `datavzrd` report directory): Entry point for the interactive HTML report.
- **`<project_id>-qc_summary.tsv`**: Project-level PASS/WARN/FAIL counts and failure rates stratified by `alignment_qc`, `variants_qc`, and `total`.
- **`<project_id>-qc_flags.tsv`**: Merged per-sample flag matrix containing `N_warning_flags`, `N_fail_flags`, overall `QC_status`, and all active flag columns.
- **`<project_id>-qc_flags_summary.tsv`**: Per-flag PASS/FAIL counts used by the report's summary bar plots.
- **`<project_id>-datavzrd.yml`**: Unified `datavzrd` report specification.
