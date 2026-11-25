# QC Criteria & Flag Rules

Cohort-level `PASS`, `WARN`, and `FAIL` classifications are driven by two YAML files:

- `--alignment_qc_criteria` (default example: [`assets/align_qc_flags.yml`](https://github.com/HTGenomeAnalysisUnit/nf-WGS-QC/blob/main/assets/align_qc_flags.yml))
- `--variant_qc_criteria` (default example: [`assets/var_qc_flags.yml`](https://github.com/HTGenomeAnalysisUnit/nf-WGS-QC/blob/main/assets/var_qc_flags.yml))

## YAML Structure

Each file defines two optional top-level mappings, `fail_rules` and `warning_rules`. Keys are rule names and values are Python-style boolean expressions referencing flattened metric column names (`<section>-<subkey>-...`):

```yaml
fail_rules:
  rule_name: "<metric_expression>"

warning_rules:
  warn_rule_name: "<metric_expression>"
```

When a rule expression evaluates to `True` for a sample, that sample is flagged (`1` in `*.qc_flags.tsv`).

### Overall Sample `QC_status`

- **`FAIL`**: One or more `fail_rules` triggered (in alignment QC, variant QC, or either when merged).
- **`WARN`**: Zero `fail_rules` triggered and one or more `warning_rules` triggered.
- **`PASS`**: Zero `fail_rules` and zero `warning_rules` triggered.

---

## Default Alignment QC Criteria (`assets/align_qc_flags.yml`)

```yaml
fail_rules:
  discordant_sex: "somalier-concordant_sex == 0"
  low_global_min_10: "mosdepth-dist_global-total-pct_bases_min_depth_10 < 0.85"
  low_region_min_10: "mosdepth-dist_region-total-pct_bases_min_depth_10 < 0.85"
  low_baseq: "samtools_stats-average_quality < 35"
  low_coverage_autosomes: "mosdepth-summary-autosomes-mean_cov < 20"
  autosome_has_zero_reads: "samtools_idxstats-autosomes-min_mapped_norm == 0"
  autosome_low_uniformity: "samtools_idxstats-autosomes-lowest_cov_fraction_of_mean < 0.7"

warning_rules:
  duplicate_samples_detected: "somalier-duplicated_samples != ''"
  unexpected_related_samples: "somalier-unexpected_related != ''"
  uncommon_karyotype: "somalier-inferred_karyotype not in ['XY', 'XX']"
```

---

## Default Variant QC Criteria (`assets/var_qc_flags.yml`)

```yaml
fail_rules:
  autosome_has_zero_vars: "bcftools_idxstats-autosomes-min_n_vars == 0"
  contaminated: "charr-charr > 0.04 or charr-inconsistent_ab_het_rate > 0.15"
  ts_tv_ratio: "bcftools_stats-tstv-ts_tv < 1.8 or bcftools_stats-tstv-ts_tv > 2.4"

warning_rules:
  var_low_coverage: "bcftools_stats-dp-fraction_genotypes_gt_10 < 0.7"
```

!!! tip "Finding exact metric names"
    For the complete catalog of available flattened metric keys for both alignment and variant QC rules, see the [QC Reports Schema](../qc_reports_schema.md).
