#!/usr/bin/env python
"""Combine variant QC manifests into a selected TSV and a detailed JSON report.

JSON stores nested reports per sample. TSV metrics join nested keys with dashes.
Large metric tables are filtered in Polars; --workers controls concurrent file
reads. Runtime dependencies are Polars, NumPy and PyYAML.
"""

import argparse
import ast
import csv
import io
import json
import math
import operator
import os
import tokenize
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from functools import reduce
from itertools import islice

import numpy as np
import polars as pl
import yaml

VERSION = 1.0

# Read QC_REPORT_PRECISION env var or set 5 as default
PRECISION = int(os.getenv("QC_REPORT_PRECISION", 5))

# Default is set as human autosomes and sex chromosomes
AUTOSOMES = [f'{i}' for i in range(1, 23)]
SEX_CHROMS = ['X', 'Y']

BCFTOOLS_STATS_BASE_KEYS = [
    'sn-number_of_snps',
    'sn-number_of_indels',
    'sn-number_of_multiallelic_sites',
    'tstv-ts',
    'tstv-tv',
    'tstv-ts_tv',
]

# Define the set of metrics to report in the summary TSV table
TSV_SELECTED_KEYS = {
    'bcftools_idxstats': {
        'visible': [
            'autosomes-n_vars',
            'autosomes-mean_n_vars_norm',
            'autosomes-min_n_vars',
            'autosomes-min_n_vars_chrom',
        ],
        'hidden': [],
    },
    'bcftools_stats': {
        'visible': list(BCFTOOLS_STATS_BASE_KEYS),
        # QUAL and DP report thresholds are added by selected_tsv_columns.
        'hidden': [],
    },
    'bcftools_stats_regions': {
        'visible': list(BCFTOOLS_STATS_BASE_KEYS),
        # QUAL and DP report thresholds are added by selected_tsv_columns.
        'hidden': [],
    },
    'charr': {
        'visible': [
            'charr',
            'heterozygosity_rate',
            'inconsistent_ab_het_rate',
        ],
        'hidden': [],
    },
}

KEY_TRANSLATION = str.maketrans('', '', '()\\/')


def log(msg: str = "", level: str = "INFO"):
    print(f"[{datetime.now().isoformat()}] - {level} - {msg}")


def sanitize_key(key: str):
    # Change spaces into underscore, remove brackets and symbols
    clean_key = key.rstrip(':').translate(KEY_TRANSLATION).replace(' ', '_').lower()
    return clean_key.replace('%', 'pct').replace('&', 'and')


def csv_list(value: str) -> list[str]:
    items = [x.strip() for x in value.split(",") if x.strip()]
    if not items:
        raise argparse.ArgumentTypeError("Expected comma-separated values")
    return items


def parse_threshold(value: str):
    try:
        numeric = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"Expected a numeric threshold, got {value!r}") from error
    if not math.isfinite(numeric) or numeric < 0:
        raise argparse.ArgumentTypeError("Threshold must be a finite non-negative number")
    return int(numeric) if numeric.is_integer() else numeric


def normalize_chrom(chrom: str) -> str:
    """Remove an exact chr prefix without stripping characters from contig names."""
    return chrom.removeprefix('chr')


def parse_value(value):
    if value is None or value in ('', 'NA', '.'):
        return None
    try:
        return int(value)
    except ValueError:
        try:
            parsed = float(value)
            return round(parsed, PRECISION) if math.isfinite(parsed) else None
        except ValueError:
            return value


def json_values(value):
    # Missing numeric measurements must be JSON null rather than NaN/Infinity.
    if isinstance(value, dict):
        return {key: json_values(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_values(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def parse_bcftools_idxstats(file_path, autosomes, canonical_chroms):
    """Parse bcftools index --stats output and compute chromosome group summaries."""
    autosomes_norm = [normalize_chrom(c) for c in autosomes]
    canonical_norm = [normalize_chrom(c) for c in canonical_chroms]
    keep = set(canonical_norm) | {'M', 'MT'}

    idxstats = {}
    total_n_vars = 0
    extra_n_vars = 0
    saw_chr_prefix = False
    saw_any_contig = False

    with open(file_path, 'r') as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3 or parts[0] == '*':
                continue
            chrom = parts[0]
            size = int(parts[1]) if parts[1] not in ('.', '', 'NA') else 0
            n_vars = int(parts[2])

            saw_any_contig = True
            if chrom.startswith('chr'):
                saw_chr_prefix = True
            total_n_vars += n_vars

            norm_chrom = normalize_chrom(chrom)
            if norm_chrom in keep:
                idxstats[chrom] = {
                    'n_vars': n_vars,
                    'n_vars_norm': round(n_vars / (size / 1e6), PRECISION) if size else 0,
                }
            else:
                extra_n_vars += n_vars

    normalized = {normalize_chrom(c): values for c, values in idxstats.items()}
    # bcftools index --stats omits contigs with zero records; populate missing canonical contigs with 0.
    for orig_chrom, norm_chrom in zip(canonical_chroms, canonical_norm):
        if norm_chrom not in normalized:
            chrom_key = (f'chr{norm_chrom}' if saw_chr_prefix else norm_chrom) if saw_any_contig else orig_chrom
            zero_entry = {'n_vars': 0, 'n_vars_norm': 0}
            idxstats[chrom_key] = zero_entry
            normalized[norm_chrom] = zero_entry

    autosome_entries = [(c, normalized[c]) for c in autosomes_norm if c in normalized]
    autosomes_n_vars = sum(v['n_vars'] for _, v in autosome_entries)
    minimum = min(autosome_entries, key=lambda item: item[1]['n_vars'], default=None)
    mean_norm = (
        round(sum(v['n_vars_norm'] for _, v in autosome_entries) / len(autosome_entries), PRECISION)
        if autosome_entries else None
    )

    idxstats['autosomes'] = {
        'n_vars': autosomes_n_vars,
        'min_n_vars': minimum[1]['n_vars'] if minimum else None,
        'min_n_vars_chrom': minimum[0] if minimum else None,
        'mean_n_vars_norm': mean_norm,
    }
    idxstats['canonical'] = {
        'n_vars': sum(normalized[c]['n_vars'] for c in canonical_norm if c in normalized),
    }
    idxstats['extra_contigs'] = {'n_vars': extra_n_vars}
    idxstats['total'] = {'n_vars': total_n_vars}
    return idxstats


def parse_charr_report(file_path, vcf_sample_id=None):
    """Parse sceVCF / CHARR TSV output, skipping leading ## comment lines."""
    with open(file_path, 'r') as handle:
        data_lines = (line for line in handle if line.strip() and not line.startswith('##'))
        rows = csv.reader(data_lines, delimiter='\t')
        header = next(rows, None)
        if header is None:
            raise ValueError(f'Empty CHARR report: {file_path}')

        cleaned_header = [col.lstrip('#').strip() for col in header]
        sample_col_idx = None
        selected_cols = []
        for idx, col in enumerate(cleaned_header):
            sanitized = sanitize_key(col)
            if sanitized == 'sample':
                sample_col_idx = idx
            else:
                selected_cols.append((idx, sanitized))

        matched_row = None
        last_row = None
        row_count = 0
        for row in rows:
            if not row:
                continue
            row_count += 1
            last_row = row
            if (
                vcf_sample_id is not None
                and sample_col_idx is not None
                and sample_col_idx < len(row)
                and row[sample_col_idx] == vcf_sample_id
            ):
                matched_row = row

    if row_count == 0:
        raise ValueError(f'No data rows in CHARR report: {file_path}')
    if matched_row is not None:
        target_row = matched_row
    elif row_count == 1:
        target_row = last_row
    else:
        raise ValueError(f'Sample {vcf_sample_id!r} not found in multi-sample CHARR report: {file_path}')

    return {key: parse_value(target_row[idx]) for idx, key in selected_cols if idx < len(target_row)}


def _parse_dp_bin(bin_str: str) -> float:
    """Convert bcftools stats DP bin labels (e.g. '10', '>500', '<0') to comparable floats."""
    bin_str = bin_str.strip()
    if bin_str.startswith('>'):
        return float(bin_str[1:]) + 0.5
    if bin_str.startswith('<'):
        return float(bin_str[1:]) - 0.5
    return float(bin_str)


def _cumulative_gt_counts(values: np.ndarray, counts: np.ndarray, thresholds) -> dict:
    """Compute sum of counts where values > threshold for each threshold."""
    if values.size == 0:
        return {t: 0 for t in thresholds}
    if values.size > 1 and np.any(values[:-1] > values[1:]):
        order = np.argsort(values, kind='mergesort')
        values = values[order]
        counts = counts[order]
    suffix = np.empty(counts.size + 1, dtype=np.int64)
    suffix[-1] = 0
    suffix[:-1] = counts[::-1].cumsum()[::-1]
    indices = np.searchsorted(values, thresholds, side='right')
    return {t: int(suffix[idx]) for t, idx in zip(thresholds, indices)}


def parse_bcftools_stats(file_path, var_qual_report=20, var_dp_report=10):
    """Parse SN, TSTV, QUAL, DP and ST sections from a sample bcftools stats file."""
    sn = {}
    tstv = {}
    st = {}

    qual_vals = []
    qual_snps = []
    qual_indels = []
    max_finite_qual = 0.0

    dp_vals = []
    dp_genos = []

    with open(file_path, 'r') as handle:
        for line in handle:
            if not line or line[0] == '#':
                continue
            if line.startswith('QUAL\t'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 7 or parts[1] != '0':
                    continue
                q_str = parts[2].strip()
                if q_str in ('.', 'NA', ''):
                    q_val = -1.0
                else:
                    q_val = float(q_str)
                    if math.isfinite(q_val) and q_val > max_finite_qual:
                        max_finite_qual = q_val
                qual_vals.append(q_val)
                qual_snps.append(int(parts[3]))
                qual_indels.append(int(parts[6]))
            elif line.startswith('DP\t'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 4 or parts[1] != '0':
                    continue
                dp_vals.append(_parse_dp_bin(parts[2]))
                dp_genos.append(int(parts[3]))
            elif line.startswith('SN\t'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 4 or parts[1] != '0':
                    continue
                sn[sanitize_key(parts[2])] = parse_value(parts[3].strip())
            elif line.startswith('TSTV\t'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 8 or parts[1] != '0':
                    continue
                ts = int(parts[2])
                tv = int(parts[3])
                ts_1st = int(parts[5])
                tv_1st = int(parts[6])
                tstv = {
                    'ts': ts,
                    'tv': tv,
                    'ts_tv': round(ts / tv, PRECISION) if tv else parse_value(parts[4].strip()),
                    'ts_1st_alt': ts_1st,
                    'tv_1st_alt': tv_1st,
                    'ts_tv_1st_alt': round(ts_1st / tv_1st, PRECISION) if tv_1st else parse_value(parts[7].strip()),
                }
            elif line.startswith('ST\t'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 4 or parts[1] != '0':
                    continue
                st[parts[2].strip()] = int(parts[3])

    # QUAL > x using x bins from 0 to max QUAL by 10 (plus configured var_qual_report)
    qual_arr = np.asarray(qual_vals, dtype=np.float64)
    snps_arr = np.asarray(qual_snps, dtype=np.int64)
    indels_arr = np.asarray(qual_indels, dtype=np.int64)
    total_qual_snps = int(snps_arr.sum())
    total_qual_indels = int(indels_arr.sum())

    qual_bins = list(range(0, int(max_finite_qual) + 1, 10))
    qual_thresholds = sorted(set(qual_bins) | {var_qual_report})
    snps_gt = _cumulative_gt_counts(qual_arr, snps_arr, qual_thresholds)
    indels_gt = _cumulative_gt_counts(qual_arr, indels_arr, qual_thresholds)

    qual = {}
    for x in qual_thresholds:
        n_snp = snps_gt[x]
        n_indel = indels_gt[x]
        qual[f'n_snps_gt_{x}'] = n_snp
        qual[f'fraction_snps_gt_{x}'] = round(n_snp / total_qual_snps, PRECISION) if total_qual_snps else None
        qual[f'n_indels_gt_{x}'] = n_indel
        qual[f'fraction_indels_gt_{x}'] = round(n_indel / total_qual_indels, PRECISION) if total_qual_indels else None

    # DP > x using x bins from 0 to 100 by 10 and > 100 all together (plus configured var_dp_report)
    dp_arr = np.asarray(dp_vals, dtype=np.float64)
    genos_arr = np.asarray(dp_genos, dtype=np.int64)
    total_genotypes = int(genos_arr.sum())

    dp_thresholds = sorted(set(range(0, 101, 10)) | {var_dp_report})
    genos_gt = _cumulative_gt_counts(dp_arr, genos_arr, dp_thresholds)

    dp = {}
    for x in dp_thresholds:
        n_geno = genos_gt[x]
        dp[f'n_genotypes_gt_{x}'] = n_geno
        dp[f'fraction_genotypes_gt_{x}'] = round(n_geno / total_genotypes, PRECISION) if total_genotypes else None

    return {
        'sn': sn,
        'tstv': tstv,
        'qual': qual,
        'dp': dp,
        'st': st,
    }


def qc_column_index(columns):
    """Tokenize the schema once into a trie for longest column-name matching."""
    aliases, trie = {}, {}
    for index, column in enumerate(columns):
        alias = f'__qc_column_{index}'
        aliases[alias] = column
        node = trie
        for token in tokenize.generate_tokens(io.StringIO(column).readline):
            if token.type not in (tokenize.NEWLINE, tokenize.ENDMARKER):
                node = node.setdefault(token.string, {})
        node[None] = alias
    return aliases, trie


def qc_condition_expr(condition, columns, *, column_index=None):
    """Translate Python conditions to Polars, accepting flattened column names.

    Known column names (including hyphenated names) take precedence over
    arithmetic. Strings are left intact. Calls and attribute access are rejected.
    """
    aliases, trie = column_index if column_index is not None else qc_column_index(columns)

    tokens = list(tokenize.generate_tokens(io.StringIO(condition).readline))
    rewritten = []
    index = 0
    while index < len(tokens):
        node, end, match = trie, index, None
        while end < len(tokens) and tokens[end].string in node:
            node = node[tokens[end].string]
            end += 1
            if None in node:
                match = (end, node[None])
        if match is not None:
            index, alias = match
            rewritten.append((tokenize.NAME, alias))
        else:
            rewritten.append((tokens[index].type, tokens[index].string))
            index += 1

    binary_ops = {
        ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
        ast.Div: operator.truediv, ast.FloorDiv: operator.floordiv,
        ast.Mod: operator.mod, ast.Pow: operator.pow,
        ast.BitAnd: operator.and_, ast.BitOr: operator.or_,
    }
    comparisons = {
        ast.Lt: operator.lt, ast.LtE: operator.le, ast.Gt: operator.gt,
        ast.GtE: operator.ge, ast.Eq: operator.eq, ast.NotEq: operator.ne,
    }

    def is_missing_marker(node):
        return isinstance(node, ast.Constant) and node.value in ('', 'NA')

    def convert(node):
        if isinstance(node, ast.Name):
            if node.id not in aliases:
                raise ValueError(f'Unknown QC column: {node.id}')
            return pl.col(aliases[node.id])
        if isinstance(node, ast.Constant):
            return pl.lit(node.value)
        if isinstance(node, ast.BoolOp):
            operation = operator.and_ if isinstance(node.op, ast.And) else operator.or_
            return reduce(operation, (convert(value) for value in node.values))
        if isinstance(node, ast.UnaryOp):
            operations = {ast.Not: operator.invert, ast.Invert: operator.invert,
                          ast.USub: operator.neg, ast.UAdd: operator.pos}
            if type(node.op) in operations:
                return operations[type(node.op)](convert(node.operand))
        if isinstance(node, ast.BinOp) and type(node.op) in binary_ops:
            return binary_ops[type(node.op)](convert(node.left), convert(node.right))
        if isinstance(node, ast.Compare):
            results = []
            left_node = node.left
            for operation, right_node in zip(node.ops, node.comparators):
                left = convert(left_node)
                if isinstance(operation, (ast.In, ast.NotIn)):
                    values = ast.literal_eval(right_node)
                    if not isinstance(values, (list, tuple, set)):
                        raise ValueError('Membership requires a literal list, tuple or set')
                    result = left.is_in(list(values))
                    if isinstance(operation, ast.NotIn):
                        result = ~result
                elif isinstance(operation, (ast.Is, ast.IsNot)):
                    if not isinstance(right_node, ast.Constant) or right_node.value is not None:
                        raise ValueError('Only "is None" and "is not None" are supported')
                    result = left.is_null()
                    if isinstance(operation, ast.IsNot):
                        result = ~result
                elif (isinstance(operation, (ast.Eq, ast.NotEq))
                      and (is_missing_marker(right_node) != is_missing_marker(left_node))):
                    # col == 'NA' / col == '' (or != ) test for a missing value: null, empty or 'NA'
                    other = left if is_missing_marker(right_node) else convert(right_node)
                    text = other.cast(pl.String, strict=False).str.strip_chars()
                    result = other.is_null() | text.is_in(['', 'NA'])
                    if isinstance(operation, ast.NotEq):
                        result = ~result
                elif type(operation) in comparisons:
                    result = comparisons[type(operation)](left, convert(right_node))
                else:
                    raise ValueError(f'Unsupported comparison: {type(operation).__name__}')
                results.append(result)
                left_node = right_node
            return reduce(operator.and_, results)
        raise ValueError(f'Unsupported QC expression: {type(node).__name__}')

    return convert(ast.parse(tokenize.untokenize(rewritten), mode='eval').body)


def evaluate_qc_flags(df, qc_flags):
    """Evaluate each rule once; null results do not trigger a flag."""
    rules = []
    conditions = []
    column_index = qc_column_index(df.columns)
    for severity in ('fail_rules', 'warning_rules'):
        for flag, condition in (qc_flags.get(severity) or {}).items():
            name = f'__qc_rule_{len(rules)}'
            try:
                expression = qc_condition_expr(condition, df.columns, column_index=column_index)
            except (ValueError, TypeError, SyntaxError, tokenize.TokenError, pl.exceptions.PolarsError) as error:
                raise ValueError(f'Invalid {severity} flag {flag!r} ({condition!r}): {error}') from error
            conditions.append(expression.alias(name))
            rules.append((severity, flag, name))

    if not rules:
        return df.with_columns(
            pl.lit('').alias('QC_fail_flags'),
            pl.lit('').alias('QC_warn_flags'),
            pl.lit(0, dtype=pl.UInt32).alias('N_FAIL_conditions'),
            pl.lit(0, dtype=pl.UInt32).alias('N_WARN_conditions'),
            pl.lit('PASS').alias('QC_status'),
        )

    try:
        evaluated = df.with_columns(conditions).select([name for _, _, name in rules])
        for rule, dtype in zip(rules, evaluated.dtypes):
            if dtype != pl.Boolean:
                raise ValueError(f'QC condition {rule[1]!r} must return booleans')
        evaluated = evaluated.fill_null(False)
    except (ValueError, TypeError, pl.exceptions.PolarsError) as error:
        raise ValueError(f'Could not evaluate QC conditions: {error}') from error

    def flag_string(severity):
        expressions = [
            pl.when(pl.col(name)).then(pl.lit(flag)).otherwise(None)
            for group, flag, name in rules if group == severity
        ]
        return pl.concat_str(expressions, separator='|', ignore_nulls=True).fill_null('') if expressions else pl.lit('')

    flags = evaluated.select(
        flag_string('fail_rules').alias('QC_fail_flags'),
        flag_string('warning_rules').alias('QC_warn_flags'),
        *[
            pl.sum_horizontal([
                pl.col(name).cast(pl.UInt32) for group, _, name in rules if group == severity
            ] or [pl.lit(0, dtype=pl.UInt32)]).alias(column)
            for severity, column in (('fail_rules', 'N_FAIL_conditions'),
                                     ('warning_rules', 'N_WARN_conditions'))
        ],
    ).with_columns(
        pl.when(pl.col('N_FAIL_conditions') > 0).then(pl.lit('FAIL'))
        .when(pl.col('N_WARN_conditions') > 0).then(pl.lit('WARN'))
        .otherwise(pl.lit('PASS')).alias('QC_status')
    )
    return df.with_columns(flags.get_columns())


def read_table(file_path, columns=None):
    # Keep identifiers such as 001 as strings in manifests and reports.
    return pl.read_csv(file_path, separator='\t', infer_schema=False, columns=columns)


def parse_metric_file(task, args):
    """Parse one file without modifying shared state (also used by workers)."""
    report, sample_id, vcf_sample_id, path = task
    canonical = args.autosomes + args.sex_chroms
    try:
        if report == 'bcftools_idxstats':
            values = parse_bcftools_idxstats(path, args.autosomes, canonical)
        elif report in ('bcftools_stats', 'bcftools_stats_regions'):
            values = parse_bcftools_stats(path, args.var_qual_report, args.var_dp_report)
        elif report == 'charr':
            values = parse_charr_report(path, vcf_sample_id=vcf_sample_id)
        else:
            raise ValueError(f'Unknown report type: {report}')
    except (OSError, ValueError, KeyError, pl.exceptions.PolarsError) as error:
        raise ValueError(f'Could not parse {report} for sample {sample_id!r} ({path}): {error}') from error
    return report, sample_id, vcf_sample_id, values


def collect_metrics(args):
    """Bound pending file reads and merge results in manifest order."""
    manifests = {
        'bcftools_idxstats': args.input_bcftools_idxstats,
        'bcftools_stats': args.input_bcftools_stats,
    }
    if args.input_bcftools_stats_regions:
        manifests['bcftools_stats_regions'] = args.input_bcftools_stats_regions
    manifests['charr'] = args.input_charr

    columns = [args.sample_id_col, args.vcf_sample_id_col, 'output_file']
    metrics = {}
    processed = 0

    def tasks():
        for report, path in manifests.items():
            manifest = read_table(path, columns).select(columns)
            log(f'Reading {manifest.height} {report} metric files')
            for sample_id, vcf_id, input_file in manifest.iter_rows():
                if args.verbose:
                    log(f'Processing {report} for sample {sample_id}, input file {input_file}')
                yield report, sample_id, vcf_id, input_file

    def merge(results):
        nonlocal processed
        for report, sample_id, vcf_id, values in results:
            result = metrics.setdefault(sample_id, {'vcf_sample_id': vcf_id})
            result[report] = values
            processed += 1
            if processed % 1000 == 0:
                log(f'Parsed {processed} metric files')

    task_iterator = tasks()
    if args.workers == 1:
        merge(parse_metric_file(task, args) for task in task_iterator)
    else:
        # Executor.map eagerly submits its input on older Python versions.
        # Small batches bound queued tasks and out-of-order parsed results.
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            while batch := list(islice(task_iterator, args.workers * 4)):
                merge(executor.map(lambda task: parse_metric_file(task, args), batch))

    return metrics


def metric_paths(values, path=()):
    for key, value in values.items():
        parts = (*path, str(key))
        if isinstance(value, dict):
            yield from metric_paths(value, parts)
        else:
            yield '-'.join(parts), parts


def build_metrics_frame(metrics, selected_columns, qc_flags):
    """Materialize only output columns and extra columns referenced by QC rules.

    Discover the complete schema before inference: metrics first appearing late
    in a large cohort must not be silently dropped by Polars' inference window.
    """
    paths = {}
    for values in metrics.values():
        paths.update(metric_paths(values))
    column_index = qc_column_index(['sample_id', *paths])
    required = set(selected_columns)
    for severity in ('fail_rules', 'warning_rules'):
        for flag, condition in (qc_flags.get(severity) or {}).items():
            try:
                expression = qc_condition_expr(condition, paths, column_index=column_index)
            except (ValueError, TypeError, SyntaxError, tokenize.TokenError) as error:
                raise ValueError(f'Invalid {severity} flag {flag!r} ({condition!r}): {error}') from error
            required.update(expression.meta.root_names())
    selected = {column: path for column, path in paths.items() if column in required}

    def rows():
        for sample_id, values in metrics.items():
            row = {'sample_id': sample_id, 'vcf_sample_id': values.get('vcf_sample_id')}
            for column, path in selected.items():
                value = values
                for key in path:
                    if not isinstance(value, dict):
                        value = None
                        break
                    value = value.get(key)
                row[column] = None if isinstance(value, float) and not math.isfinite(value) else value
            yield row

    return pl.from_dicts(rows(), infer_schema_length=None,
                         schema_overrides={'sample_id': pl.String, 'vcf_sample_id': pl.String})


def write_metrics_json(metrics, output_path):
    """Normalize and encode one sample at a time to limit temporary memory."""
    with open(output_path, 'w') as handle:
        handle.write('{')
        for index, (sample_id, values) in enumerate(metrics.items()):
            handle.write(',\n' if index else '\n')
            handle.write('    ' + json.dumps(sample_id) + ': ')
            encoded = json.dumps(json_values(values), indent=4, allow_nan=False)
            handle.write(encoded.replace('\n', '\n    '))
        handle.write('\n}' if metrics else '}')


def parse_arguments(argv=None):
    # Handle --version early so it works without required arguments
    version_parser = argparse.ArgumentParser(add_help=False)
    version_parser.add_argument("--version", action="version", version=f"{VERSION}")
    version_parser.parse_known_args(argv)

    parser = argparse.ArgumentParser(description="Summarise variant QC")
    parser.add_argument('--workers', type=int, default=1, help='Concurrent metric file readers (default: 1; increase for I/O-bound storage)')
    parser.add_argument('--verbose', action='store_true', help='Log every metric file')
    parser.add_argument("--sample_id_col", help="Column name for sample ID in input files", default="id")
    parser.add_argument("--vcf_sample_id_col", help="Column name for VCF sample ID in input files", default="vcf_sample_id")
    parser.add_argument("--qc_flags", required=True, help="QC flags YAML file")
    parser.add_argument("--input_bcftools_idxstats", required=True, help="Input bcftools idxstats manifest file")
    parser.add_argument("--input_bcftools_stats", required=True, help="Input bcftools stats manifest file")
    parser.add_argument(
        "--input_bcftools_stats_regions",
        "--input_bcftools_stats_region",
        dest="input_bcftools_stats_regions",
        default=None,
        help="Optional input bcftools stats on regions manifest file",
    )
    parser.add_argument(
        "--input_charr",
        "--input_charr_contamination",
        dest="input_charr",
        required=True,
        help="Input CHARR contamination manifest file",
    )
    parser.add_argument("--autosomes", type=csv_list, default=AUTOSOMES, help="Comma-separated list of autosomes")
    parser.add_argument("--sex_chroms", type=csv_list, default=SEX_CHROMS, help="Comma-separated list of sex chromosomes")
    parser.add_argument("--var_qual_report", type=parse_threshold, default=20, help="QUAL threshold for reported variant counts/fractions (default: 20)")
    parser.add_argument("--var_dp_report", type=parse_threshold, default=10, help="DP threshold for reported genotype counts/fractions (default: 10)")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files (TSV, JSON, summary tables)")

    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    return args


def selected_tsv_columns(var_qual_report, var_dp_report, include_regions=False):
    columns = ['sample_id', 'vcf_sample_id']
    stats_threshold_keys = [
        f'qual-n_snps_gt_{var_qual_report}',
        f'qual-fraction_snps_gt_{var_qual_report}',
        f'qual-n_indels_gt_{var_qual_report}',
        f'qual-fraction_indels_gt_{var_qual_report}',
        f'dp-n_genotypes_gt_{var_dp_report}',
        f'dp-fraction_genotypes_gt_{var_dp_report}',
    ]
    tools = ['bcftools_idxstats', 'bcftools_stats']
    if include_regions:
        tools.append('bcftools_stats_regions')
    tools.append('charr')

    for tool in tools:
        categories = TSV_SELECTED_KEYS[tool]
        columns.extend(f'{tool}-{key}' for key in categories['visible'] + categories['hidden'])
        if tool in ('bcftools_stats', 'bcftools_stats_regions'):
            columns.extend(f'{tool}-{key}' for key in stats_threshold_keys)
    return columns


def write_tabular_reports(df, qc_flags, accepted_tsv_columns, output_prefix):
    """Write the two-header metric table, flag matrix and rule descriptions."""
    qc_tail_columns = ['QC_fail_flags', 'QC_warn_flags', 'N_FAIL_conditions', 'N_WARN_conditions']
    id_columns = [col for col in ('sample_id', 'vcf_sample_id') if col in df.columns]
    metric_columns = [
        col for col in accepted_tsv_columns
        if col in df.columns and col not in id_columns and col not in qc_tail_columns and col != 'QC_status'
    ]
    tail_columns = [col for col in qc_tail_columns if col in df.columns]
    df = df.select(id_columns + ['QC_status'] + metric_columns + tail_columns)

    # Polars does not support pandas-style MultiIndex headers, so write the two-row TSV header manually.
    output_columns = []
    for col in df.columns:
        if col in ('sample_id', 'vcf_sample_id'):
            output_columns.append((col, 'sampleid_section'))
        elif col.startswith('bcftools_stats_regions-'):
            _, suffix = col.split('-', 1)
            output_columns.append((f'regions-{suffix}', 'bcftools_stats_regions'))
        else:
            parts = col.split('-', 1)
            if len(parts) == 2:
                output_columns.append((parts[1], parts[0]))
            else:
                output_columns.append((col, 'qc_section'))

    primary_header, secondary_header = zip(*output_columns)
    output_path = f"{output_prefix}.tsv"
    with open(output_path, 'w', newline='') as f:
        f.write('\t'.join(primary_header) + '\n')
        f.write('\t'.join(secondary_header) + '\n')
    with open(output_path, 'ab') as f:
        # Treat empty/whitespace-only strings as missing so they are written as NA like nulls (df itself is left untouched)
        df_out = df.with_columns(
            pl.when(pl.col(pl.String).str.strip_chars() == "")
            .then(None)
            .otherwise(pl.col(pl.String))
            .name.keep()
        )
        df_out.write_csv(f, separator='\t', include_header=False, null_value='NA', quote_style='never')
    log(f"TSV file saved to {output_prefix}.tsv")

    # Export one indicator column per flag along with counts and overall status.
    flag_names_warn = qc_flags.get('warning_rules') or {}
    flag_names_fail = qc_flags.get('fail_rules') or {}
    qc_flags_status_df = df.select(
        'sample_id', 'vcf_sample_id',
        'QC_status',
        'N_WARN_conditions',
        'N_FAIL_conditions',
        *[pl.col('QC_warn_flags').str.split('|').list.contains(flag).cast(pl.UInt8).alias(f"WARN-{flag}")
          for flag in flag_names_warn],
        *[pl.col('QC_fail_flags').str.split('|').list.contains(flag).cast(pl.UInt8).alias(f"FAIL-{flag}")
          for flag in flag_names_fail],
    )
    qc_flags_status_df.write_csv(f"{output_prefix}.qc_flags.tsv", separator='\t')
    log(f"QC flags status file saved to {output_prefix}.qc_flags.tsv")

    # Dump the QC criteria YML as a TSV table with columns for severity, flag name and description
    qc_criteria_rows = []
    for severity, flags in qc_flags.items():
        for flag, description in (flags or {}).items():
            qc_criteria_rows.append({
                'severity': severity,
                'flag': flag,
                'description': description,
            })
    qc_criteria_df = pl.DataFrame(qc_criteria_rows, schema={
        'severity': pl.String, 'flag': pl.String, 'description': pl.String})
    qc_criteria_df.write_csv(f"{output_prefix}.qc_criteria.tsv", separator='\t')
    log(f"QC criteria file saved to {output_prefix}.qc_criteria.tsv")

    return df, qc_flags_status_df, output_columns


def build_datavzrd_config(df, qc_flags_status_df, output_columns, output_prefix):
    """Describe the variant QC report tables for datavzrd."""
    secondary_header = [section for _, section in output_columns]
    sections = list(dict.fromkeys(secondary_header))
    palette = ["#eeeeee", "#c6dbef", "#c7e9c0", "#fdd0a2", "#dadaeb", "#fcbba1", "#d9d9d9", "#bcbddc", "#fde0dd"]
    section_colors = [palette[i % len(palette)] for i in range(len(sections))]

    datavzrd_columns = {}
    for col, (primary, _) in zip(df.columns, output_columns):
        if col == 'sample_id':
            datavzrd_columns[primary] = {'label': 'Sample ID', 'display-mode': 'pinned'}
        elif col == 'vcf_sample_id':
            datavzrd_columns[primary] = {'label': 'VCF Sample ID', 'display-mode': 'pinned'}
        elif col == 'QC_status':
            datavzrd_columns[primary] = {
                'label': 'QC status',
                'display-mode': 'pinned',
                'plot': {
                    'heatmap': {
                        'scale': 'ordinal',
                        'domain': ['PASS', 'WARN', 'FAIL'],
                        'range': ['#00a000', '#ffa500', '#ff0000'],
                    }
                },
            }
        else:
            datavzrd_columns.setdefault(primary, {'label': primary, 'precision': 3})

    datavzrd_config = {
        'name': 'QC metrics report',
        'datasets': {
            'var_qc_metrics': {
                'path': f"{output_prefix}.tsv",
                'separator': '\t',
                'headers': 2,
            }
        },
        'views': {
            'var_qc_metrics': {
                'dataset': 'var_qc_metrics',
                'desc': 'Variant QC metrics',
                'render-table': {
                    'headers': {
                        1: {
                            'label': 'QC section',
                            'plot': {
                                'heatmap': {
                                    'scale': 'ordinal',
                                    'domain': sections,
                                    'range': section_colors,
                                }
                            },
                        }
                    },
                    'columns': datavzrd_columns,
                },
            }
        },
    }

    # Add the QC flags status table: flags with value 1 in red, QC_status colored by outcome
    flag_columns = [c for c in qc_flags_status_df.columns if c.startswith(('WARN-', 'FAIL-'))]
    qc_flags_columns = {
        'sample_id': {'label': 'Sample ID', 'display-mode': 'pinned'},
        'vcf_sample_id': {'label': 'VCF Sample ID', 'display-mode': 'pinned'},
        'QC_status': {
            'display-mode': 'pinned',
            'plot': {
                'heatmap': {
                    'scale': 'ordinal',
                    'domain': ['PASS', 'WARN', 'FAIL'],
                    'range': ['#00a000', '#ffa500', '#ff0000'],
                }
            },
        },
        'N_WARN_conditions': {'display-mode': 'pinned'},
        'N_FAIL_conditions': {'display-mode': 'pinned'},
    }
    for c in flag_columns:
        qc_flags_columns[c] = {
            'plot': {
                'heatmap': {
                    'scale': 'ordinal',
                    'domain': ['0', '1'],
                    'range': ['#ffffff', '#ff0000'],
                }
            }
        }
    datavzrd_config['datasets']['var_qc_flags'] = {
        'path': f"{output_prefix}.qc_flags.tsv",
        'separator': '\t',
    }
    datavzrd_config['views']['var_qc_flags'] = {
        'dataset': 'var_qc_flags',
        'desc': 'Variant QC flags status (1 = flag raised)',
        'render-table': {'columns': qc_flags_columns},
    }

    # Add the QC criteria table: severity colored (fail_rules red, warning_rules orange)
    datavzrd_config['datasets']['var_qc_criteria'] = {
        'path': f"{output_prefix}.qc_criteria.tsv",
        'separator': '\t',
    }
    datavzrd_config['views']['var_qc_criteria'] = {
        'dataset': 'var_qc_criteria',
        'desc': 'Variant QC criteria',
        'render-table': {
            'columns': {
                'severity': {
                    'label': 'Severity',
                    'plot': {
                        'heatmap': {
                            'scale': 'ordinal',
                            'domain': ['warning_rules', 'fail_rules'],
                            'range': ['#ffa500', '#ff0000'],
                        }
                    },
                },
                'flag': {'label': 'Flag'},
                'description': {'label': 'Description'},
            }
        },
    }

    return datavzrd_config


def main(argv=None):
    args = parse_arguments(argv)
    log(f'Sample ID column: {args.sample_id_col}')
    log(f'VCF Sample ID column: {args.vcf_sample_id_col}')
    log(f'Canonical chromosomes: {args.autosomes + args.sex_chroms}')
    log(f'Variant QUAL report threshold: {args.var_qual_report}')
    log(f'Variant DP report threshold: {args.var_dp_report}')
    log(f'File readers: {args.workers}')
    accepted_tsv_columns = selected_tsv_columns(
        args.var_qual_report,
        args.var_dp_report,
        include_regions=bool(args.input_bcftools_stats_regions),
    )

    with open(args.qc_flags, 'r') as f:
        qc_flags = yaml.safe_load(f) or {}
    log(f"QC flags: {qc_flags}")

    qc_metrics = collect_metrics(args)
    if not qc_metrics:
        raise ValueError('No samples found in input manifests')
    log(f'Finished reading metrics for {len(qc_metrics)} samples')
    write_metrics_json(qc_metrics, f'{args.output_prefix}.json')
    log(f'JSON file saved to {args.output_prefix}.json')
    df = build_metrics_frame(qc_metrics, accepted_tsv_columns, qc_flags)
    del qc_metrics

    df = evaluate_qc_flags(df, qc_flags)

    df, flags_df, output_columns = write_tabular_reports(df, qc_flags, accepted_tsv_columns, args.output_prefix)
    datavzrd_config = build_datavzrd_config(df, flags_df, output_columns, args.output_prefix)
    with open(f"{args.output_prefix}.qc_var_datavzrd.yml", 'w') as f:
        yaml.safe_dump(datavzrd_config, f, sort_keys=False)
    log(f"datavzrd config saved to {args.output_prefix}.qc_var_datavzrd.yml")


if __name__ == '__main__':
    main()
