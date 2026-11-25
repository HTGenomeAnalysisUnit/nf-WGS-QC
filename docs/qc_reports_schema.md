# QC Summary Reports Schema (`summarise_align_qc.py` & `summarise_var_qc.py`)

This document describes the output file schemas produced by:

- `bin/summarise_align_qc.py` — Alignment-level QC summary
- `bin/summarise_var_qc.py` — Variant-level QC summary

---

## 1. Shared Conventions

### Generated Output Files

Given `--output_prefix <PREFIX>`, both scripts generate five files:

| File | Format | Description |
| :--- | :--- | :--- |
| `<PREFIX>.json` | JSON | Nested per-sample dictionary containing all collected detailed metrics. |
| `<PREFIX>.tsv` | 2-header TSV | Selected summary metrics per sample, annotated with `QC_status` and rule flags. |
| `<PREFIX>.qc_flags.tsv` | TSV | Binary indicator matrix (`0`/`1`) per sample for every configured `WARN-*` and `FAIL-*` rule. |
| `<PREFIX>.qc_criteria.tsv` | TSV | Table of configured QC rules (`severity`, `flag`, `description`). |
| `<PREFIX>.qc_align_datavzrd.yml` / `<PREFIX>.qc_var_datavzrd.yml` | YAML | `datavzrd` report configuration for rendering the three TSV tables above. |

### Formatting, Precision & Missing Values

- **Sample keys**: Top-level keys in `<PREFIX>.json` (and `sample_id` in TSV files) correspond to the user sample ID (`--sample_id_col`, default `id`). Identifiers are always preserved as strings (e.g. `"001"`).
- **Floating-point precision**: Controlled by the `QC_REPORT_PRECISION` environment variable (default: `5` decimal places).
- **Missing / non-finite values**: Empty strings, `.`, `NA`, `NaN`, and `±Inf` are encoded as `null` in JSON (`allow_nan=False`) and `NA` in `<PREFIX>.tsv`.

### Two-Header TSV Layout & Flattened QC Rule Names

`<PREFIX>.tsv` uses a **two-row header** designed for `datavzrd`:

1. **Row 1 (Primary Header)**: Column identifier within its section (used as the unique column key in `datavzrd`).
2. **Row 2 (Secondary Header)**: QC section name (`sampleid_section`, `qc_section`, or the tool/manifest section name), rendered as a color-coded header bar in `datavzrd`.

When writing filter expressions in QC YAML files (e.g. `assets/align_qc_flags.yml`, `assets/var_qc_flags.yml`), reference metrics using their **flattened dash-joined JSON path** (`<section>-<subkey>-...`):

```yaml
fail_rules:
  contaminated: "charr-charr > 0.05"
  low_dp10: "bcftools_stats-dp-fraction_genotypes_gt_10 < 0.85"

warning_rules:
  low_tstv: "bcftools_stats-tstv-ts_tv < 2.0"
```

- **Any JSON metric can be used in QC rules**: Both scripts use the same `build_metrics_frame` + `qc_condition_expr` engine, which automatically materializes any metric path present in `<PREFIX>.json` (including special names such as `samtools_stats-1st_fragments` or `bcftools_stats-st-C>T`) when referenced by `fail_rules` or `warning_rules`, even if that metric is not included in `<PREFIX>.tsv`. Note that a referenced path must exist in `<PREFIX>.json` for **at least one sample** in the run; referencing an optional section omitted for the entire cohort (e.g. `bcftools_stats_regions-*` without `--input_bcftools_stats_regions`) raises an `Unknown QC column` error.
- **Null handling & operators in rules**:
  - Polars `null` (`null` in JSON) propagates through arithmetic, comparisons (`<`, `<=`, `>`, `>=`, `==`, `!=`), and membership (`in`, `not in`) and is converted to `False` (`fill_null(False)`), so missing measurements never trigger a flag unless explicitly tested—even if the column is `null` across all samples.
  - Explicit missingness checks: `col == 'NA'` and `col == ''` match `null`, empty/whitespace strings (`""`, e.g. empty `somalier-duplicated_samples`), and `'NA'` (with `!= 'NA'` / `!= ''` negating this); `col is None` and `col is not None` test strictly for `null`.
  - Supported operators: arithmetic (`+`, `-`, `*`, `/`, `//`, `%`, `**`, unary `+`/`-`), comparisons (`<`, `<=`, `>`, `>=`, `==`, `!=`), membership (`in`, `not in` against literal lists/tuples/sets), and boolean logic (`and`, `or`, `not`, or `&`, `|`, `~` — note that `&`, `|`, `~` follow Python bitwise precedence and require parentheses around comparisons, e.g. `(a > 1) | (b > 2)`).
- **Status precedence**:
  - `FAIL` if `N_FAIL_conditions > 0`
  - `WARN` if `N_FAIL_conditions == 0` and `N_WARN_conditions > 0`
  - `PASS` otherwise

### Companion Tables (`*.qc_flags.tsv` & `*.qc_criteria.tsv`)

- **`<PREFIX>.qc_flags.tsv`** (single header row):
  - `sample_id` (`string`), `bam_sample_id` or `vcf_sample_id` (`string`)
  - `QC_status` (`string`: `PASS` | `WARN` | `FAIL`)
  - `N_WARN_conditions` (`uint32`), `N_FAIL_conditions` (`uint32`)
  - `WARN-<flag>` (`uint8`: `0` or `1`) for each rule in `warning_rules`
  - `FAIL-<flag>` (`uint8`: `0` or `1`) for each rule in `fail_rules`
- **`<PREFIX>.qc_criteria.tsv`** (single header row):
  - `severity` (`string`: `fail_rules` | `warning_rules`), `flag` (`string`), `description` (`string`: rule expression)

---

## 2. Alignment QC (`summarise_align_qc.py`)

### 2.1 JSON Output Schema (`*-align_qc.json`)

Top-level object keyed by `sample_id`:

```json
{
  "<sample_id>": {
    "bam_sample_id": "<string>",
    "samtools_stats": { ... },
    "samtools_idxstats": { ... },
    "mosdepth": {
      "summary": { ... },
      "dist_global": { ... },
      "dist_region": { ... }
    },
    "verifybamid": { ... },
    "qrs": { ... },
    "somalier": { ... }
  }
}
```

#### Section Details

- **`bam_sample_id`** (`string`): Sample identifier in the BAM header (`--bam_sample_id_col`, default `bam_sample_id`).
- **`samtools_stats`** (`object`):
  - All `SN` summary numbers from `samtools stats`, key-sanitized to lowercase snake_case (`float`), e.g. `raw_total_sequences`, `filtered_sequences`, `sequences`, `is_sorted`, `1st_fragments`, `last_fragments`, `reads_mapped`, `reads_mapped_and_paired`, `reads_unmapped`, `reads_properly_paired`, `reads_paired`, `reads_duplicated`, `reads_mq0`, `reads_qc_failed`, `non-primary_alignments`, `supplementary_alignments`, `total_length`, `bases_mapped`, `bases_mapped_cigar`, `bases_duplicated`, `mismatches`, `error_rate`, `average_length`, `maximum_length`, `average_quality`, `insert_size_average`, `insert_size_standard_deviation`, `inward_oriented_pairs`, `outward_oriented_pairs`, `pairs_with_other_orientation`, `pairs_on_different_chromosomes`, `percentage_of_properly_paired_reads_pct`.
  - Derived ratios (`float | null`, `[0, 1]`):
    - `reads_pct_dup_on_all`: `reads_duplicated / sequences`
    - `reads_pct_dup_on_mapped`: `reads_duplicated / reads_mapped`
    - `reads_pct_mapped`: `reads_mapped / sequences`
    - `reads_pct_properly_paired`: `reads_properly_paired / sequences`
  - Base quality thresholds for each `Q` in `--target_base_qualities` (default `30`, `float | null`, `[0, 1]`):
    - `bases_at_least_Q<Q>_FFQ`: fraction of first-fragment bases with base quality `>= Q`
    - `bases_at_least_Q<Q>_LFQ`: fraction of last-fragment bases with base quality `>= Q`
    - `bases_at_least_Q<Q>`: combined fraction of bases with base quality `>= Q`
- **`samtools_idxstats`** (`object`):
  - Per-contig objects for canonical chromosomes (`--autosomes` + `--sex_chroms`) and mitochondrial (`M`/`MT` if present), e.g. `"chr1"`, ..., `"chrX"`, `"chrY"`, `"chrM"`:
    - `mapped` (`int`), `mapped_norm` (`float`, mapped reads per Mb), `unmapped` (`int`), `fraction_mapped` (`float`)
    - Sex chromosomes (`--sex_chroms`) additionally include `scaled_cov` (`float | null`): `mapped_norm` scaled against the first autosome (`mapped_norm / autosome1_mapped_norm * 2`).
  - `autosomes` (`object`):
    - `mapped` (`int`), `unmapped` (`int`), `fraction_mapped` (`float`)
    - `min_mapped_norm` (`float | null`): lowest `mapped_norm` across autosomes
    - `min_mapped_norm_chrom` (`string | null`): normalized autosome name (without `chr`) with lowest `mapped_norm`
    - `lowest_cov_fraction_of_mean` (`float`): `min_mapped_norm / mean_autosome_mapped_norm`
  - `canonical` (`object`):
    - `mapped` (`int`), `unmapped` (`int`), `fraction_mapped` (`float`)
- **`mosdepth`** (`object`):
  - `summary` (`object`):
    - Per-contig entries (`chr1`, ..., `chrM`, `total`, plus `*_region` and `total_region` when region BED was used): `bases` (`int`), `mean` (`float`), `min` (`int`), `max` (`int`).
    - Group aggregates (`autosomes`, `canonical`, `extra_contigs`, plus `autosomes_region`, `canonical_region` when region BED was used): `bases` (`int`), `mean_cov` (`float | null`, unweighted mean across contigs), `min_cov` (`int | null`), `max_cov` (`int | null`).
  - `dist_global` and (when `.mosdepth.region.dist.txt` is provided) `dist_region` (`object`):
    - Per-contig (`chr1`, ..., `total`) and group (`autosomes`, `canonical`) objects containing `pct_bases_min_depth_<D>` (`float`, `[0, 1]`) for each depth `D` in `--target_depths` (default `10,20,30`).
- **`verifybamid`** (`object`):
  - `freemix` (`float | null`): DNA contamination estimate from VerifyBamID2 (`FREEMIX`).
- **`qrs`** (`object`):
  - From `*.per_sample_run.tsv` (when present): `n_runs` (`int`), `n_flowcells` (`int`), `n_lanes` (`int`).
  - From `*.read_stats.tsv` (when present): `mean_baseq` (`float`), `median_baseq` (`float`), `mean_read_gc_perc` (`float`).
- **`somalier`** (`object`):
  - Ancestry: `ancestry` (`string`), `ancestry_prob` (`float`).
  - Sex & karyotype: `provided_sex` (`string`), `X_depth_mean` (`float`), `Y_depth_mean` (`float`), `scaled_X` (`float | null`), `scaled_Y` (`float | null`), `X_n_vars` (`int`), `X_het` (`int`), `Y_n_vars` (`int`), `inferred_sex` (`string`: `male` | `female` | `unknown`), `inferred_karyotype` (`string`: `XX` | `XY` | `XXY` | `XXX` | `X` | `unknown`), `concordant_sex` (`int`: `1` if concordant or `provided_sex == 'unknown'`, else `0`).
  - Relatedness: `duplicated_samples` (`string`: `|`-delimited sample IDs with relatedness `>= --relatedness_duplicate_value`, empty `""` if none), `unexpected_related` (`string`: `|`-delimited sample IDs with unexpected relatedness, empty `""` if none).

### 2.2 Main TSV Output Schema (`*-align_qc.tsv`)

Columns below assume defaults `--target_base_qualities 30` and `--target_depths 10,20,30`. Optional columns depend on whether the corresponding input reports (e.g. `*.mosdepth.region.dist.txt` or `*.read_stats.tsv`) are present in the manifests.

| Row 1 (Primary Header) | Row 2 (Secondary Header) | Flattened Name (for `align_qc_flags.yml`) | Type | Description |
| :--- | :--- | :--- | :--- | :--- |
| `sample_id` | `sampleid_section` | `sample_id` | `string` | User sample ID |
| `bam_sample_id` | `sampleid_section` | `bam_sample_id` | `string` | Sample ID in BAM header |
| `QC_status` | `qc_section` | — *(output)* | `string` | `PASS`, `WARN`, or `FAIL` |
| `sequences` | `samtools_stats` | `samtools_stats-sequences` | `float` | Total read sequences |
| `reads_mapped` | `samtools_stats` | `samtools_stats-reads_mapped` | `float` | Mapped reads |
| `reads_properly_paired` | `samtools_stats` | `samtools_stats-reads_properly_paired` | `float` | Properly paired reads |
| `reads_duplicated` | `samtools_stats` | `samtools_stats-reads_duplicated` | `float` | Duplicated reads |
| `reads_mq0` | `samtools_stats` | `samtools_stats-reads_mq0` | `float` | Reads with mapping quality 0 |
| `total_length` | `samtools_stats` | `samtools_stats-total_length` | `float` | Total sequenced bases |
| `error_rate` | `samtools_stats` | `samtools_stats-error_rate` | `float` | Mismatches per mapped base |
| `average_length` | `samtools_stats` | `samtools_stats-average_length` | `float` | Mean read length |
| `average_quality` | `samtools_stats` | `samtools_stats-average_quality` | `float` | Mean base quality |
| `reads_pct_dup_on_all` | `samtools_stats` | `samtools_stats-reads_pct_dup_on_all` | `float` | Fraction duplicated over all reads |
| `reads_pct_dup_on_mapped` | `samtools_stats` | `samtools_stats-reads_pct_dup_on_mapped` | `float` | Fraction duplicated over mapped reads |
| `reads_pct_mapped` | `samtools_stats` | `samtools_stats-reads_pct_mapped` | `float` | Fraction mapped reads |
| `reads_pct_properly_paired` | `samtools_stats` | `samtools_stats-reads_pct_properly_paired` | `float` | Fraction properly paired reads |
| `bases_at_least_Q<Q>_FFQ` | `samtools_stats` | `samtools_stats-bases_at_least_Q<Q>_FFQ` | `float` | Fraction R1 bases with `QUAL >= Q` |
| `bases_at_least_Q<Q>_LFQ` | `samtools_stats` | `samtools_stats-bases_at_least_Q<Q>_LFQ` | `float` | Fraction R2 bases with `QUAL >= Q` |
| `bases_at_least_Q<Q>` | `samtools_stats` | `samtools_stats-bases_at_least_Q<Q>` | `float` | Fraction total bases with `QUAL >= Q` |
| `autosomes-mapped` | `samtools_idxstats` | `samtools_idxstats-autosomes-mapped` | `int` | Total mapped reads on autosomes |
| `autosomes-min_mapped_norm` | `samtools_idxstats` | `samtools_idxstats-autosomes-min_mapped_norm` | `float` | Lowest mapped reads/Mb across autosomes |
| `autosomes-min_mapped_norm_chrom` | `samtools_idxstats` | `samtools_idxstats-autosomes-min_mapped_norm_chrom` | `string` | Autosome with lowest mapped reads/Mb |
| `autosomes-lowest_cov_fraction_of_mean` | `samtools_idxstats` | `samtools_idxstats-autosomes-lowest_cov_fraction_of_mean` | `float` | Lowest autosome norm mapped / autosome mean |
| `autosomes-fraction_mapped` | `samtools_idxstats` | `samtools_idxstats-autosomes-fraction_mapped` | `float` | Mapped / (mapped + unmapped) on autosomes |
| `canonical-mapped` | `samtools_idxstats` | `samtools_idxstats-canonical-mapped` | `int` | Total mapped reads on canonical chroms |
| `canonical-fraction_mapped` | `samtools_idxstats` | `samtools_idxstats-canonical-fraction_mapped` | `float` | Mapped / (mapped + unmapped) on canonical chroms |
| `summary-autosomes-bases` | `mosdepth` | `mosdepth-summary-autosomes-bases` | `int` | Sequenced bases on autosomes |
| `summary-autosomes-mean_cov` | `mosdepth` | `mosdepth-summary-autosomes-mean_cov` | `float` | Mean coverage across autosomes |
| `summary-canonical-bases` | `mosdepth` | `mosdepth-summary-canonical-bases` | `int` | Sequenced bases on canonical chroms |
| `summary-canonical-mean_cov` | `mosdepth` | `mosdepth-summary-canonical-mean_cov` | `float` | Mean coverage across canonical chroms |
| `summary-extra_contigs-bases` | `mosdepth` | `mosdepth-summary-extra_contigs-bases` | `int` | Sequenced bases on extra contigs |
| `summary-extra_contigs-mean_cov` | `mosdepth` | `mosdepth-summary-extra_contigs-mean_cov` | `float` | Mean coverage across extra contigs |
| `summary-autosomes_region-mean_cov`\* | `mosdepth` | `mosdepth-summary-autosomes_region-mean_cov` | `float` | Mean coverage on target regions (autosomes) |
| `summary-canonical_region-mean_cov`\* | `mosdepth` | `mosdepth-summary-canonical_region-mean_cov` | `float` | Mean coverage on target regions (canonical) |
| `dist_region-total-pct_bases_min_depth_<D>`\* | `mosdepth` | `mosdepth-dist_region-total-pct_bases_min_depth_<D>` | `float` | Fraction target-region bases at depth `>= D` |
| `dist_global-total-pct_bases_min_depth_<D>` | `mosdepth` | `mosdepth-dist_global-total-pct_bases_min_depth_<D>` | `float` | Fraction genome bases at depth `>= D` |
| `freemix` | `verifybamid` | `verifybamid-freemix` | `float` | VerifyBamID2 contamination estimate |
| `median_baseq`\* | `qrs` | `qrs-median_baseq` | `float` | Median read base quality (`*.read_stats.tsv`) |
| `mean_read_gc_perc`\* | `qrs` | `qrs-mean_read_gc_perc` | `float` | Mean read GC fraction (`*.read_stats.tsv`) |
| `n_runs`\* | `qrs` | `qrs-n_runs` | `int` | Number of sequencing runs (`*.per_sample_run.tsv`) |
| `n_flowcells`\* | `qrs` | `qrs-n_flowcells` | `int` | Number of flowcells (`*.per_sample_run.tsv`) |
| `n_lanes`\* | `qrs` | `qrs-n_lanes` | `int` | Number of lanes (`*.per_sample_run.tsv`) |
| `ancestry` | `somalier` | `somalier-ancestry` | `string` | Predicted ancestry label |
| `ancestry_prob` | `somalier` | `somalier-ancestry_prob` | `float` | Probability of predicted ancestry |
| `provided_sex` | `somalier` | `somalier-provided_sex` | `string` | Pedigree/manifest sex (`male`, `female`, `unknown`) |
| `scaled_X` | `somalier` | `somalier-scaled_X` | `float` | Normalized chrX depth (`X_depth_mean / depth_mean * 2`) |
| `scaled_Y` | `somalier` | `somalier-scaled_Y` | `float` | Normalized chrY depth (`Y_depth_mean / depth_mean * 2`) |
| `X_n_vars` | `somalier` | `somalier-X_n_vars` | `int` | Assayed chrX variant sites |
| `X_het` | `somalier` | `somalier-X_het` | `int` | Heterozygous chrX calls |
| `Y_n_vars` | `somalier` | `somalier-Y_n_vars` | `int` | Assayed chrY variant sites |
| `inferred_sex` | `somalier` | `somalier-inferred_sex` | `string` | Inferred sex (`male`, `female`, `unknown`) |
| `inferred_karyotype` | `somalier` | `somalier-inferred_karyotype` | `string` | Inferred sex karyotype (`XX`, `XY`, `XXY`, `XXX`, `X`, `unknown`) |
| `concordant_sex` | `somalier` | `somalier-concordant_sex` | `int` | `1` if `provided_sex` matches `inferred_sex` (or `unknown`), else `0` |
| `duplicated_samples` | `somalier` | `somalier-duplicated_samples` | `string` | `|`-separated duplicate sample IDs (`NA` in TSV if none) |
| `unexpected_related` | `somalier` | `somalier-unexpected_related` | `string` | `|`-separated unexpected related sample IDs (`NA` in TSV if none) |
| `QC_fail_flags` | `qc_section` | — *(output)* | `string` | `|`-separated triggered `fail_rules` flags (`NA` if none) |
| `QC_warn_flags` | `qc_section` | — *(output)* | `string` | `|`-separated triggered `warning_rules` flags (`NA` if none) |
| `N_FAIL_conditions` | `qc_section` | — *(output)* | `uint32` | Number of triggered `fail_rules` |
| `N_WARN_conditions` | `qc_section` | — *(output)* | `uint32` | Number of triggered `warning_rules` |

*\*Included when the corresponding optional report files are present in the input manifests.*

---

## 3. Variant QC (`summarise_var_qc.py`)

### 3.1 JSON Output Schema (`*-var_qc.json`)

Top-level object keyed by `sample_id`:

```json
{
  "<sample_id>": {
    "vcf_sample_id": "<string>",
    "bcftools_idxstats": {
      "chr1": { "n_vars": 363238, "n_vars_norm": 1459.0425 },
      "...": {},
      "autosomes": {
        "n_vars": 4635765,
        "min_n_vars": 71004,
        "min_n_vars_chrom": "22",
        "mean_n_vars_norm": 1607.69506
      },
      "canonical": { "n_vars": 4758493 },
      "extra_contigs": { "n_vars": 0 },
      "total": { "n_vars": 4758506 }
    },
    "bcftools_stats": {
      "sn": { ... },
      "tstv": { ... },
      "qual": { ... },
      "dp": { ... },
      "st": { ... }
    },
    "bcftools_stats_regions": {
      "sn": { ... },
      "tstv": { ... },
      "qual": { ... },
      "dp": { ... },
      "st": { ... }
    },
    "charr": { ... }
  }
}
```

#### Section Details

- **`vcf_sample_id`** (`string`): Sample identifier in the VCF file (`--vcf_sample_id_col`, default `vcf_sample_id`).
- **`bcftools_idxstats`** (`object`):
  - Per-chromosome entries for every canonical chromosome (`--autosomes` + `--sex_chroms`) and mitochondrial (`M`/`MT` if present), e.g. `"chr1"`, ..., `"chr22"`, `"chrX"`, `"chrY"`, `"chrM"`:
    - `n_vars` (`int`): Variant count on the chromosome (canonical chromosomes absent from `bcftools index --stats` default to `0`).
    - `n_vars_norm` (`float`): Variants per Mb (`n_vars / (chrom_length / 1e6)`).
  - `autosomes` (`object`):
    - `n_vars` (`int`): Total variants across all autosomes.
    - `min_n_vars` (`int | null`): Variant count in the autosome with the lowest `n_vars`.
    - `min_n_vars_chrom` (`string | null`): Normalized name (without `chr`) of the autosome with the lowest `n_vars`.
    - `mean_n_vars_norm` (`float | null`): Unweighted mean of `n_vars_norm` across autosomes.
  - `canonical` (`object`): `n_vars` (`int`) — Total variants across autosomes + sex chromosomes.
  - `extra_contigs` (`object`): `n_vars` (`int`) — Total variants on non-canonical, non-mitochondrial contigs.
  - `total` (`object`): `n_vars` (`int`) — Total variants across all contigs in the index.
- **`bcftools_stats`** and *(optional, when `--input_bcftools_stats_regions` is provided)* **`bcftools_stats_regions`** (`object`):
  - `sn` (`object`): All `SN` summary values for set `0` (`int`), sanitized to lowercase snake_case:
    - `number_of_samples`, `number_of_records`, `number_of_no-alts`, `number_of_snps`, `number_of_mnps`, `number_of_indels`, `number_of_others`, `number_of_multiallelic_sites`, `number_of_multiallelic_snp_sites`.
  - `tstv` (`object`): Transitions/transversions for set `0`:
    - `ts` (`int`), `tv` (`int`), `ts_tv` (`float | null`), `ts_1st_alt` (`int`), `tv_1st_alt` (`int`), `ts_tv_1st_alt` (`float | null`).
  - `qual` (`object`): Cumulative counts and fractions strictly above threshold `x` (`QUAL > x`) for `x` in `0, 10, 20, ..., max_qual_bin` (plus `--var_qual_report`, default `20`). Denominators are total SNPs and indels in the `QUAL` distribution:
    - `n_snps_gt_<x>` (`int`), `fraction_snps_gt_<x>` (`float | null`, `[0, 1]`)
    - `n_indels_gt_<x>` (`int`), `fraction_indels_gt_<x>` (`float | null`, `[0, 1]`)
  - `dp` (`object`): Cumulative genotype counts and fractions strictly above depth `x` (`DP > x`) for `x` in `0, 10, 20, ..., 100` (where `x = 100` aggregates all genotypes with `DP > 100`, including `>500`, plus `--var_qual_dp`, default `10`):
    - `n_genotypes_gt_<x>` (`int`), `fraction_genotypes_gt_<x>` (`float | null`, `[0, 1]`)
  - `st` (`object`): All 12 nucleotide substitution counts (`int`):
    - `A>C`, `A>G`, `A>T`, `C>A`, `C>G`, `C>T`, `G>A`, `G>C`, `G>T`, `T>A`, `T>C`, `T>G`.
- **`charr`** (`object`):
  - All per-sample CHARR (`sceVCF`) columns normalized to lowercase (`int | float | null`):
    - `hom_total`, `hom_covered`, `hq_hom`, `hq_hom_rate`, `het_total`, `het_covered`, `hq_het`, `hq_het_rate`, `charr`, `mean_ref_ab_hom_alt`, `heterozygosity_rate`, `inconsistent_ab_het_rate`.

### 3.2 Main TSV Output Schema (`*-var_qc.tsv`)

Columns appear in the exact deterministic order below (shown with defaults `--var_qual_report 20` and `--var_qual_dp 10`; substituting `<Q>` and `<D>` if overridden). For `bcftools_stats_regions`, Row 1 primary headers are prefixed with `regions-` so all `datavzrd` column keys remain unique.

| # | Row 1 (Primary Header) | Row 2 (Secondary Header) | Flattened Name (for `var_qc_flags.yml`) | Type | Description |
| :- | :--- | :--- | :--- | :--- | :--- |
| 1 | `sample_id` | `sampleid_section` | `sample_id` | `string` | User sample ID |
| 2 | `vcf_sample_id` | `sampleid_section` | `vcf_sample_id` | `string` | Sample ID in VCF file |
| 3 | `QC_status` | `qc_section` | — *(output)* | `string` | `PASS`, `WARN`, or `FAIL` |
| 4 | `autosomes-n_vars` | `bcftools_idxstats` | `bcftools_idxstats-autosomes-n_vars` | `int` | Total variants in autosomes |
| 5 | `autosomes-mean_n_vars_norm` | `bcftools_idxstats` | `bcftools_idxstats-autosomes-mean_n_vars_norm` | `float` | Mean variants/Mb across autosomes |
| 6 | `autosomes-min_n_vars` | `bcftools_idxstats` | `bcftools_idxstats-autosomes-min_n_vars` | `int` | Variant count in autosome with lowest count |
| 7 | `autosomes-min_n_vars_chrom` | `bcftools_idxstats` | `bcftools_idxstats-autosomes-min_n_vars_chrom` | `string` | Autosome name with lowest variant count |
| 8 | `sn-number_of_snps` | `bcftools_stats` | `bcftools_stats-sn-number_of_snps` | `int` | Number of SNPs |
| 9 | `sn-number_of_indels` | `bcftools_stats` | `bcftools_stats-sn-number_of_indels` | `int` | Number of indels |
| 10 | `sn-number_of_multiallelic_sites` | `bcftools_stats` | `bcftools_stats-sn-number_of_multiallelic_sites` | `int` | Number of multiallelic sites |
| 11 | `tstv-ts` | `bcftools_stats` | `bcftools_stats-tstv-ts` | `int` | Number of transitions |
| 12 | `tstv-tv` | `bcftools_stats` | `bcftools_stats-tstv-tv` | `int` | Number of transversions |
| 13 | `tstv-ts_tv` | `bcftools_stats` | `bcftools_stats-tstv-ts_tv` | `float` | Transitions / transversions ratio (`ts / tv`) |
| 14 | `qual-n_snps_gt_20` | `bcftools_stats` | `bcftools_stats-qual-n_snps_gt_20` | `int` | SNPs with `QUAL > --var_qual_report` |
| 15 | `qual-fraction_snps_gt_20` | `bcftools_stats` | `bcftools_stats-qual-fraction_snps_gt_20` | `float` | Fraction of SNPs with `QUAL > --var_qual_report` |
| 16 | `qual-n_indels_gt_20` | `bcftools_stats` | `bcftools_stats-qual-n_indels_gt_20` | `int` | Indels with `QUAL > --var_qual_report` |
| 17 | `qual-fraction_indels_gt_20` | `bcftools_stats` | `bcftools_stats-qual-fraction_indels_gt_20` | `float` | Fraction of indels with `QUAL > --var_qual_report` |
| 18 | `dp-n_genotypes_gt_10` | `bcftools_stats` | `bcftools_stats-dp-n_genotypes_gt_10` | `int` | Genotypes with `DP > --var_qual_dp` |
| 19 | `dp-fraction_genotypes_gt_10` | `bcftools_stats` | `bcftools_stats-dp-fraction_genotypes_gt_10` | `float` | Fraction of genotypes with `DP > --var_qual_dp` |
| 20\* | `regions-sn-number_of_snps` | `bcftools_stats_regions` | `bcftools_stats_regions-sn-number_of_snps` | `int` | Number of SNPs in target regions |
| 21\* | `regions-sn-number_of_indels` | `bcftools_stats_regions` | `bcftools_stats_regions-sn-number_of_indels` | `int` | Number of indels in target regions |
| 22\* | `regions-sn-number_of_multiallelic_sites` | `bcftools_stats_regions` | `bcftools_stats_regions-sn-number_of_multiallelic_sites` | `int` | Number of multiallelic sites in target regions |
| 23\* | `regions-tstv-ts` | `bcftools_stats_regions` | `bcftools_stats_regions-tstv-ts` | `int` | Number of transitions in target regions |
| 24\* | `regions-tstv-tv` | `bcftools_stats_regions` | `bcftools_stats_regions-tstv-tv` | `int` | Number of transversions in target regions |
| 25\* | `regions-tstv-ts_tv` | `bcftools_stats_regions` | `bcftools_stats_regions-tstv-ts_tv` | `float` | `ts / tv` ratio in target regions |
| 26\* | `regions-qual-n_snps_gt_20` | `bcftools_stats_regions` | `bcftools_stats_regions-qual-n_snps_gt_20` | `int` | Target-region SNPs with `QUAL > --var_qual_report` |
| 27\* | `regions-qual-fraction_snps_gt_20` | `bcftools_stats_regions` | `bcftools_stats_regions-qual-fraction_snps_gt_20` | `float` | Fraction target-region SNPs with `QUAL > --var_qual_report` |
| 28\* | `regions-qual-n_indels_gt_20` | `bcftools_stats_regions` | `bcftools_stats_regions-qual-n_indels_gt_20` | `int` | Target-region indels with `QUAL > --var_qual_report` |
| 29\* | `regions-qual-fraction_indels_gt_20` | `bcftools_stats_regions` | `bcftools_stats_regions-qual-fraction_indels_gt_20` | `float` | Fraction target-region indels with `QUAL > --var_qual_report` |
| 30\* | `regions-dp-n_genotypes_gt_10` | `bcftools_stats_regions` | `bcftools_stats_regions-dp-n_genotypes_gt_10` | `int` | Target-region genotypes with `DP > --var_qual_dp` |
| 31\* | `regions-dp-fraction_genotypes_gt_10` | `bcftools_stats_regions` | `bcftools_stats_regions-dp-fraction_genotypes_gt_10` | `float` | Fraction target-region genotypes with `DP > --var_qual_dp` |
| 32 | `charr` | `charr` | `charr-charr` | `float` | CHARR DNA contamination estimator |
| 33 | `heterozygosity_rate` | `charr` | `charr-heterozygosity_rate` | `float` | Heterozygosity rate (`HET / (HOM + HET)`) |
| 34 | `inconsistent_ab_het_rate` | `charr` | `charr-inconsistent_ab_het_rate` | `float` | Fraction of het calls with inconsistent allele balance |
| 35 | `QC_fail_flags` | `qc_section` | — *(output)* | `string` | `|`-separated triggered `fail_rules` flags (`NA` if none) |
| 36 | `QC_warn_flags` | `qc_section` | — *(output)* | `string` | `|`-separated triggered `warning_rules` flags (`NA` if none) |
| 37 | `N_FAIL_conditions` | `qc_section` | — *(output)* | `uint32` | Number of triggered `fail_rules` |
| 38 | `N_WARN_conditions` | `qc_section` | — *(output)* | `uint32` | Number of triggered `warning_rules` |

*\*Columns 20–31 are present only when `--input_bcftools_stats_regions` is provided.*
