#!/usr/bin/env python
"""Combine alignment QC manifests into a selected TSV and a detailed JSON report.

JSON stores nested reports per sample. TSV metrics join nested keys with dashes.
Large metric tables are filtered in Polars; --workers controls concurrent file
reads. Runtime dependencies are Polars, NumPy and PyYAML.
"""

import ast
import argparse
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

# Define the set of metrics to report in the table
TSV_SELECTED_KEYS = {
    'samtools_stats' : {
        'visible': [
            'sequences',
            'total_length',
            'average_length',
            'average_quality',
            'error_rate',
            'reads_pct_dup_on_all',
            'reads_pct_dup_on_mapped',
            'reads_pct_mapped',
            'reads_pct_properly_paired'
            # bases_at_least_Q{t} will be added downstream according to arguments
        ],
        'hidden': [
            'reads_mq0',
            'reads_duplicated',
            'reads_mapped',
            'reads_properly_paired'
            # bases_at_least_Q{t}_FFQ will be added downstream according to arguments
            # bases_at_least_Q{t}_LFQ will be added downstream according to arguments
        ]
    },
    'samtools_idxstats' : {
        'visible': [
            'autosomes-fraction_mapped',
            'autosomes-lowest_cov_fraction_of_mean',
            'autosomes-min_mapped_norm',
            'autosomes-min_mapped_norm_chrom',
            'canonical-fraction_mapped',
        ],
        'hidden': [
            'autosomes-mapped',
            'canonical-mapped',
        ]
    },
    'mosdepth' : {
        'visible': [
            'summary-autosomes-mean_cov',
            'summary-autosomes_region-mean_cov',
            'summary-canonical-mean_cov',
            'summary-canonical_region-mean_cov',
            'summary-total-mean_cov',
            'summary-total_region-mean_cov',
            # Distribution thresholds are added by selected_tsv_columns.
        ],
        'hidden': [
            'summary-extra_contigs-mean_cov',
            'summary-autosomes-bases',
            'summary-canonical-bases',
            'summary-extra_contigs-bases',
        ]
    },
    'somalier': {
        'visible': [
            'duplicated_samples',
            'unexpected_related',
            'provided_sex',
            'inferred_sex',
            'concordant_sex',
            'inferred_karyotype',
            'ancestry',
        ],
        'hidden': [
            'ancestry_prob',
            'scaled_X',
            'scaled_Y',
            'X_n_vars',
            'Y_n_vars',
            'X_het',
        ]
    },
    'verifybamid' : {
        'visible': ['freemix'],
        'hidden': []
    },
    'qrs' : {
        'visible': ['median_baseq', 'mean_read_gc_perc'],
        'hidden': ['n_runs', 'n_flowcells', 'n_lanes']
    }
}

KEY_TRANSLATION = str.maketrans('', '', '()\\/')


def log(msg: str = "", level: str = "INFO"):
    print(f"[{datetime.now().isoformat()}] - {level} - {msg}")


def sanitize_key(key: str):
    # Change spaces into underscore, remove brackets and symbols
    clean_key = key.rstrip(':').translate(KEY_TRANSLATION).replace(' ','_').lower()
    return clean_key.replace('%', 'pct').replace('&', 'and')


def csv_list(value: str) -> list[str]:
    items = [x.strip() for x in value.split(",") if x.strip()]
    if not items:
        raise argparse.ArgumentTypeError("Expected comma-separated values")
    return items


def parse_samtools_stats(file_path, base_q_thresholds=(30,)):
    """Read summary numbers and aggregate quality histograms across cycles."""
    base_q_thresholds = tuple(dict.fromkeys(base_q_thresholds))
    baseq_block_done = False
    stats = {}
    histograms = {section: np.zeros(0, dtype=np.int64) for section in ('FFQ', 'LFQ')}
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t', 3)
            # Parse summary stats (SN) block
            if parts[0] == 'SN':
                key = sanitize_key(parts[1])
                value = float(parts[2].strip())
                stats[key] = value
            # Parse base quality blocks (FFQ and LFQ)
            elif parts[0] in ['FFQ', 'LFQ']:
                # parts[1] is the cycle; parts[2:] contains counts for Q0, Q1, ...
                section = parts[0]
                counts = np.fromstring(line.split('\t', 2)[2], dtype=np.int64, sep='\t')
                if counts.size > histograms[section].size:
                    histograms[section] = np.pad(histograms[section], (0, counts.size - histograms[section].size))
                histograms[section][:counts.size] += counts
                baseq_block_done = True
            # Comments separate FFQ and LFQ; stop at the next data block.
            elif parts[0] and baseq_block_done:
                break
    total_bases, passing_q_bases = {}, {}
    for section, counts in histograms.items():
        total_bases[section] = int(counts.sum())
        suffix = counts[::-1].cumsum()[::-1]
        passing_q_bases[section] = {
            t: int(suffix[t]) if t < suffix.size else 0 for t in base_q_thresholds
        }
    stats['reads_pct_dup_on_all'] = stats['reads_duplicated'] / stats['sequences'] if stats['sequences'] else None
    stats['reads_pct_dup_on_mapped'] = stats['reads_duplicated'] / stats['reads_mapped'] if stats['reads_mapped'] else None
    stats['reads_pct_mapped'] = stats['reads_mapped'] / stats['sequences'] if stats['sequences'] else None
    stats['reads_pct_properly_paired'] = stats['reads_properly_paired'] / stats['sequences'] if stats['sequences'] else None
    for section in ("FFQ", "LFQ"):
        for t in base_q_thresholds:
            fraction = passing_q_bases[section][t] / total_bases[section] if total_bases[section] else float("nan")
            stats[f"bases_at_least_Q{t}_{section}"] = fraction

    n_total_bases = sum(total_bases.values())
    for t in base_q_thresholds:
        fraction = (passing_q_bases['FFQ'][t] + passing_q_bases['LFQ'][t]) / n_total_bases if n_total_bases else float("nan")
        stats[f"bases_at_least_Q{t}"] = fraction

    # Round measurements once at the parser boundary.
    for key in stats:
        if isinstance(stats[key], float):
            stats[key] = round(stats[key], PRECISION)
    return stats


def normalize_chrom(chrom):
    """Remove an exact chr prefix without stripping characters from contig names."""
    return chrom.removeprefix('chr')


def parse_samtools_idxstats(file_path, autosomes, sex_chroms, canonical_chroms):
    autosomes_norm = [normalize_chrom(c) for c in autosomes]
    canonical_norm = [normalize_chrom(c) for c in canonical_chroms]
    keep = set(canonical_norm) | {'M', 'MT'}
    idxstats = {}
    with open(file_path) as handle:
        for line in handle:
            parts = line.split('\t')
            if len(parts) < 4 or normalize_chrom(parts[0]) not in keep:
                continue
            chrom, size, mapped, unmapped = parts[:4]
            size, mapped, unmapped = int(size), int(mapped), int(unmapped)
            idxstats[chrom] = {
                'mapped': mapped,
                'mapped_norm': round(mapped / (size / 1e6) if size else 0, PRECISION),
                'unmapped': unmapped,
                'fraction_mapped': round(mapped / (mapped + unmapped) if mapped + unmapped else 0, PRECISION),
            }

    normalized = {normalize_chrom(c): values for c, values in idxstats.items()}
    for group, chroms in (('autosomes', autosomes_norm), ('canonical', canonical_norm)):
        entries = [(c, normalized[c]) for c in chroms if c in normalized]
        mapped = sum(v['mapped'] for _, v in entries)
        unmapped = sum(v['unmapped'] for _, v in entries)
        result = {'mapped': mapped, 'unmapped': unmapped}
        if group == 'autosomes':
            minimum = min(entries, key=lambda item: item[1]['mapped_norm'], default=None)
            mean = sum(v['mapped_norm'] for _, v in entries) / len(entries) if entries else 0
            result.update({
                'min_mapped_norm': minimum[1]['mapped_norm'] if minimum else None,
                'min_mapped_norm_chrom': minimum[0] if minimum else None,
                'lowest_cov_fraction_of_mean': round(minimum[1]['mapped_norm'] / mean, PRECISION) if mean else 0,
            })
        result['fraction_mapped'] = round(mapped / (mapped + unmapped), PRECISION) if mapped + unmapped else 0
        idxstats[group] = result

    scaling = normalized.get(autosomes_norm[0], {}).get('mapped_norm', 0)
    for chrom in sex_chroms:
        values = normalized.get(normalize_chrom(chrom))
        if values is not None:
            values['scaled_cov'] = round(values['mapped_norm'] / scaling * 2, PRECISION) if scaling else None
    return idxstats


def parse_simple_report(file_path, columns=()):
    """Read the final row of a single-sample report, projecting columns once."""
    with open(file_path) as handle:
        rows = csv.reader(handle, delimiter='\t')
        header = next(rows, None)
        if header is None:
            raise ValueError(f'Empty report: {file_path}')
        selected = [(i, sanitize_key(key)) for i, key in enumerate(header)
                    if not columns or key in columns]
        values = None
        for row in rows:
            if row:
                values = row
    if values is None:
        raise ValueError(f'No data rows in report: {file_path}')
    return {key: parse_value(values[i]) for i, key in selected if i < len(values)}


def aggregate_coverage(entries):
    """Preserve the report's unweighted mean across contigs."""
    return {
        'bases': sum(v['bases'] for v in entries),
        'mean_cov': round(sum(v['mean'] for v in entries) / len(entries), PRECISION) if entries else None,
        'min_cov': min((v['min'] for v in entries), default=None),
        'max_cov': max((v['max'] for v in entries), default=None),
    }


def parse_mosdepth_summary(file_path, autosomes, canonical_chroms):
    coverage = {}
    with open(file_path) as handle:
        next(handle)
        for line in handle:
            chrom, _, bases, mean, minimum, maximum = line.rstrip().split('\t')
            coverage[chrom] = {'bases': int(bases), 'mean': float(mean),
                               'min': int(minimum), 'max': int(maximum)}
    normalized = {normalize_chrom(c): values for c, values in coverage.items()}
    groups = {
        'autosomes': [normalize_chrom(c) for c in autosomes],
        'canonical': [normalize_chrom(c) for c in canonical_chroms],
    }
    keep = set(groups['canonical']) | {'M', 'MT', 'total'}
    keep.update(f'{c}_region' for c in tuple(keep))
    extras = [c for c in normalized if c not in keep]
    for group, chroms in groups.items():
        coverage[group] = aggregate_coverage([normalized[c] for c in chroms if c in normalized])
    coverage['extra_contigs'] = aggregate_coverage([normalized[c] for c in extras])

    if not any(f'{c}_region' in normalized for c in groups['canonical']):
        return coverage
    for group, chroms in groups.items():
        coverage[f'{group}_region'] = aggregate_coverage(
            [normalized[f'{c}_region'] for c in chroms if f'{c}_region' in normalized])
    for chrom in extras:
        coverage.pop(chrom, None)
        coverage.pop(f'chr{chrom}', None)
    return coverage


def parse_mosdepth_distribution(file_path, target_depths, autosomes, canonical_chroms):
    """Parse large distributions in Polars, materializing only requested metrics."""
    autosomes_norm = [normalize_chrom(c) for c in autosomes]
    canonical_norm = [normalize_chrom(c) for c in canonical_chroms]
    keep = set(canonical_norm) | {'M', 'MT', 'total'}
    keep.update(f'{c}_region' for c in tuple(keep))
    selected = (
        pl.scan_csv(file_path, separator='\t', has_header=False,
                    schema={'chrom': pl.String, 'depth': pl.Int64, 'fraction': pl.Float64})
        .filter(pl.col('depth').is_in(target_depths)
                & pl.col('chrom').str.strip_prefix('chr').is_in(keep))
        .collect()
    )
    distribution = {}
    for chrom, depth, fraction in selected.iter_rows():
        distribution.setdefault(chrom, {})[f'pct_bases_min_depth_{depth}'] = fraction
    normalized = {normalize_chrom(c): values for c, values in distribution.items()}
    for group, chroms in (('autosomes', autosomes_norm), ('canonical', canonical_norm)):
        for depth in target_depths:
            key = f'pct_bases_min_depth_{depth}'
            values = [normalized[c][key] for c in chroms if c in normalized and key in normalized[c]]
            if values:
                distribution.setdefault(group, {})[key] = round(sum(values) / len(values), PRECISION)
    return distribution


def parse_somalier_pairs(file_path, relatedness_duplicate_value, relatedness_max_discrepancy, id_map=None):
    """Filter pairs in Polars; retain only IDs and flagged relationships in Python.

    Pair tables grow quadratically with cohort size. Two projected scans avoid
    loading the full pair table and its unrelated measurements into memory.
    """
    id_map = id_map or {}
    pairs = pl.scan_csv(file_path, separator='\t', schema_overrides={
        '#sample_a': pl.String, 'sample_b': pl.String,
        'relatedness': pl.Float64, 'expected_relatedness': pl.Float64,
    })
    sample_ids = pairs.select(
        pl.concat_list('#sample_a', 'sample_b').explode().unique(maintain_order=True).alias('sample_id')
    ).collect(engine='streaming').get_column('sample_id')
    related = {id_map.get(sample, sample): {'duplicated_samples': [], 'unexpected_related': []}
               for sample in sample_ids}
    duplicate = pl.col('relatedness') >= relatedness_duplicate_value
    unexpected = (
        ((pl.col('expected_relatedness') == -1) & (pl.col('relatedness') >= 0.05))
        | ((pl.col('expected_relatedness') != -1)
           & ((pl.col('relatedness') - pl.col('expected_relatedness')).abs() >= relatedness_max_discrepancy))
    )
    flagged = pairs.select('#sample_a', 'sample_b',
                           duplicate.fill_null(False).alias('duplicate'),
                           unexpected.fill_null(False).alias('unexpected')).filter(
                               pl.col('duplicate') | pl.col('unexpected')).collect(engine='streaming')
    for sample_a, sample_b, is_duplicate, is_unexpected in flagged.iter_rows():
        sample_a, sample_b = id_map.get(sample_a, sample_a), id_map.get(sample_b, sample_b)
        for flag, active in (('duplicated_samples', is_duplicate), ('unexpected_related', is_unexpected)):
            if active:
                related[sample_a][flag].append(sample_b)
                related[sample_b][flag].append(sample_a)
    for values in related.values():
        for flag, samples in values.items():
            values[flag] = '|'.join(dict.fromkeys(samples))
    return related


def parse_somalier_samples(file_path, id_map=None):
    id_map = id_map or {}
    columns = ['sample_id', 'original_pedigree_sex', 'X_depth_mean', 'Y_depth_mean',
               'depth_mean', 'X_n', 'X_het', 'Y_n']
    samples_df = pl.read_csv(file_path, separator='\t', columns=columns,
                             schema_overrides={'sample_id': pl.String, 'original_pedigree_sex': pl.String})
    samples_df = samples_df.with_columns(*[
        pl.when(pl.col('depth_mean') > 0)
        .then(pl.col(f'{chrom}_depth_mean') / pl.col('depth_mean') * 2)
        .otherwise(None).round(PRECISION).alias(f'scaled_{chrom}')
        for chrom in ('X', 'Y')
    ])
    samples = {}
    for row in samples_df.iter_rows(named=True):
        sample_id = id_map.get(row['sample_id'], row['sample_id'])
        scaled_x, scaled_y = row['scaled_X'], row['scaled_Y']
        inferred_sex = 'unknown' if scaled_x is None or scaled_y is None else (
            'male' if (0.8 <= scaled_y <= 1.5 or (row['X_het'] < 5 and scaled_x <= 1.2)) else 'female')
        provided_sex = row['original_pedigree_sex'] or 'unknown'
        samples[sample_id] = {
            'provided_sex': provided_sex,
            'X_depth_mean': row['X_depth_mean'], 'Y_depth_mean': row['Y_depth_mean'],
            'scaled_X': scaled_x, 'scaled_Y': scaled_y,
            'X_n_vars': row['X_n'], 'X_het': row['X_het'], 'Y_n_vars': row['Y_n'],
            'inferred_sex': inferred_sex,
            'inferred_karyotype': infer_karyotype(scaled_x, scaled_y, row['X_het'], row['Y_n']),
            'concordant_sex': int(provided_sex.lower() in (inferred_sex, 'unknown')),
        }
    return samples


def parse_somalier_ancestry(file_path, id_map=None):
    id_map = id_map or {}
    ancestry_df = pl.scan_csv(file_path, separator='\t', schema_overrides={
        '#sample_id': pl.String, 'given_ancestry': pl.String,
    }).filter(pl.col('given_ancestry').is_null() | (pl.col('given_ancestry') == '')).collect()
    ancestry = {}
    for row in ancestry_df.iter_rows(named=True):
        sample_id = id_map.get(row['#sample_id'], row['#sample_id'])
        ancestry[sample_id] = {
            'ancestry': row['predicted_ancestry'],
            'ancestry_prob': row[f"{row['predicted_ancestry']}_prob"],
        }
    return ancestry

def infer_karyotype(scaled_X, scaled_Y, X_het, Y_n):
    if scaled_X is None or scaled_Y is None:
        return "unknown"
    if (1.8 <= scaled_X <= 2.2) and (scaled_Y < 0.2):
        return "XX"
    if (0.8 <= scaled_X <= 1.2) and (0.8 <= scaled_Y <= 1.5):
        return "XY"
    if (1.8 <= scaled_X <= 2.2) and (0.8 <= scaled_Y <= 1.5):
        return "XXY"
    if (2.8 <= scaled_X <= 3.2) and (scaled_Y < 0.5):
        return "XXX"
    if (0.8 <= scaled_X <= 1.2) and (Y_n < 5 and scaled_Y < 0.5):
        return "X"
    return "unknown"

def read_table(file_path, columns=None):
    # Keep identifiers such as 001 as strings in manifests and reports.
    return pl.read_csv(file_path, separator='\t', infer_schema=False, columns=columns)

def parse_value(value):
    if value is None or value in ('', 'NA', '.'):
        return None
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
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

def parse_metric_file(task, args):
    """Parse one file without modifying shared state (also used by workers)."""
    report, sample_id, bam_sample_id, path = task
    canonical = args.autosomes + args.sex_chroms
    section = None
    try:
        if report == 'samtools_stats':
            values = parse_samtools_stats(path, args.target_base_qualities)
        elif report == 'samtools_idxstats':
            values = parse_samtools_idxstats(path, args.autosomes, args.sex_chroms, canonical)
        elif report == 'mosdepth':
            if path.endswith('.mosdepth.summary.txt'):
                section = 'summary'
                values = parse_mosdepth_summary(path, args.autosomes, canonical)
            elif path.endswith(('.mosdepth.global.dist.txt', '.mosdepth.region.dist.txt')):
                section = 'dist_global' if path.endswith('.global.dist.txt') else 'dist_region'
                values = parse_mosdepth_distribution(path, args.target_depths, args.autosomes, canonical)
            else:
                raise ValueError('Unrecognised mosdepth report')
        elif report == 'verifybamid':
            values = parse_simple_report(path, columns=('FREEMIX',))
        elif report == 'qrs':
            if path.endswith('.per_sample_run.tsv'):
                values = parse_simple_report(path, columns=('N_RUNS', 'N_FLOWCELLS', 'N_LANES'))
            elif path.endswith('.read_stats.tsv'):
                values = parse_simple_report(path, columns=('MEAN_BASEQ', 'MEDIAN_BASEQ', 'MEAN_READ_GC_PERC'))
            else:
                log(f'Unrecognised QRS report: {path}', level='WARN')
                values = {}
        else:
            raise ValueError(f'Unknown report type: {report}')
    except (OSError, ValueError, KeyError, pl.exceptions.PolarsError) as error:
        raise ValueError(f'Could not parse {report} for sample {sample_id!r} ({path}): {error}') from error
    return report, sample_id, bam_sample_id, section, values


def collect_metrics(args):
    """Bound pending file reads and merge results in manifest order."""
    manifests = {
        'samtools_stats': args.input_samtools_stats,
        'samtools_idxstats': args.input_samtools_idxstats,
        'mosdepth': args.input_mosdepth,
        'verifybamid': args.input_verifybamid,
        'qrs': args.input_qrs,
    }
    columns = [args.sample_id_col, args.bam_sample_id_col, 'output_file']
    metrics, id_map = {}, {}
    processed = 0

    def tasks():
        for report, path in manifests.items():
            manifest = read_table(path, columns).select(columns)
            log(f'Reading {manifest.height} {report} metric files')
            for sample_id, bam_id, input_file in manifest.iter_rows():
                if args.verbose:
                    log(f'Processing {report} for sample {sample_id}, input file {input_file}')
                yield report, sample_id, bam_id, input_file

    def merge(results):
        nonlocal processed
        for report, sample_id, bam_id, section, values in results:
            result = metrics.setdefault(sample_id, {'bam_sample_id': bam_id})
            id_map[bam_id] = sample_id
            if report in ('samtools_stats', 'samtools_idxstats', 'verifybamid'):
                result[report] = values
            else:
                destination = result.setdefault(report, {})
                if section is not None:
                    destination = destination.setdefault(section, {})
                destination.update(values)
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

    log('Reading Somalier files')
    for results in (
        parse_somalier_ancestry(args.somalier_ancestry, id_map),
        parse_somalier_samples(args.somalier_samples, id_map),
        parse_somalier_pairs(args.somalier_pairs, args.relatedness_duplicate_value,
                             args.relatedness_max_discrepancy, id_map),
    ):
        for sample_id, values in results.items():
            metrics.setdefault(sample_id, {}).setdefault('somalier', {}).update(values)
    # A single-sample cohort has a header-only pair table. Its samples still
    # need empty relationship fields for the usual duplicate/relatedness rules.
    for values in metrics.values():
        if 'somalier' in values:
            values['somalier'].setdefault('duplicated_samples', '')
            values['somalier'].setdefault('unexpected_related', '')
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
            row = {'sample_id': sample_id, 'bam_sample_id': values.get('bam_sample_id')}
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
                         schema_overrides={'sample_id': pl.String, 'bam_sample_id': pl.String})


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
    version_parser.parse_known_args()
    
    parser = argparse.ArgumentParser(description="Summarise alignment QC")
    parser.add_argument('--workers', type=int, default=1, help='Concurrent metric file readers (default: 1; increase for I/O-bound storage)')
    parser.add_argument('--verbose', action='store_true', help='Log every metric file')
    parser.add_argument("--sample_id_col", help="Column name for sample ID in input files", default="id")
    parser.add_argument("--bam_sample_id_col", help="Column name for BAM sample ID in input files", default="bam_sample_id")
    parser.add_argument("--qc_flags", required=True, help="QC flags YAML file")
    parser.add_argument("--somalier_ancestry", required=True, help="Somalier ancestry TSV report")
    parser.add_argument("--somalier_pairs", required=True, help="Somalier pairs TSV report")
    parser.add_argument("--relatedness_duplicate_value", type=float, default=0.85, help="Pairs with relatedness equal or greater to this value are considered duplicates")
    parser.add_argument("--relatedness_max_discrepancy", type=float, default=0.2, help="Maximum allowed discrepancy between computed and expected relatedness")
    parser.add_argument("--somalier_samples", required=True, help="Somalier samples TSV report")
    parser.add_argument("--input_mosdepth", required=True, help="Input mosdepth manifest file")
    parser.add_argument("--input_samtools_idxstats", required=True, help="Input samtools idxstats manifest file")
    parser.add_argument("--input_samtools_stats", required=True, help="Input samtools stats manifest file")
    parser.add_argument("--input_verifybamid", required=True, help="Input verifybamid manifest file")
    parser.add_argument("--input_qrs", required=True, help="Input QRS manifest file")
    parser.add_argument("--autosomes", type=csv_list, default=AUTOSOMES, help="Comma-separated list of autosomes")
    parser.add_argument("--sex_chroms", type=csv_list, default=SEX_CHROMS, help="Comma-separated list of sex chromosomes")
    parser.add_argument("--target_depths", type=csv_list, default=['10', '20', '30'], help="Comma-separated list of target depths to report (default: 10,20,30)")
    parser.add_argument("--target_base_qualities", type=csv_list, default=['30'], help="Comma-separated list of target bases qualities to report (default: 30)")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files (TSV, JSON, summary tables)")
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    try:
        target_depths = [int(value) for value in args.target_depths]
        target_base_qualities = [int(value) for value in args.target_base_qualities]
    except ValueError:
        parser.error('Target depths and base qualities must be comma-separated integers')
    if any(value < 0 for value in target_depths + target_base_qualities):
        parser.error('Target depths and base qualities must be non-negative')
    args.target_depths = target_depths
    args.target_base_qualities = target_base_qualities

    return args


def selected_tsv_columns(target_base_qualities, target_depths):
    columns = ['sample_id', 'bam_sample_id']
    for tool, categories in TSV_SELECTED_KEYS.items():
        columns.extend(f'{tool}-{key}' for key in categories['visible'] + categories['hidden'])
    for quality in target_base_qualities:
        columns.extend(f'samtools_stats-bases_at_least_Q{quality}{suffix}'
                       for suffix in ('', '_FFQ', '_LFQ'))
    for depth in target_depths:
        columns.extend(f'mosdepth-{section}-total-pct_bases_min_depth_{depth}'
                       for section in ('dist_global', 'dist_region'))
    return columns


def main(argv=None):
    args = parse_arguments(argv)
    log(f'Sample ID column: {args.sample_id_col}')
    log(f'BAM Sample ID column: {args.bam_sample_id_col}')
    log(f'Canonical chromosomes: {args.autosomes + args.sex_chroms}')
    log(f'File readers: {args.workers}')
    accepted_tsv_columns = selected_tsv_columns(args.target_base_qualities, args.target_depths)

    # Read QC flags from YAML file and log them
    with open(args.qc_flags, 'r') as f:
        qc_flags = yaml.safe_load(f) or {}
    log(f"QC flags: {qc_flags}")

    qc_metrics = collect_metrics(args)
    if not qc_metrics:
        raise ValueError('No samples found in input manifests or Somalier reports')
    log(f'Finished reading metrics for {len(qc_metrics)} samples')
    write_metrics_json(qc_metrics, f'{args.output_prefix}.json')
    log(f'JSON file saved to {args.output_prefix}.json')
    df = build_metrics_frame(qc_metrics, accepted_tsv_columns, qc_flags)
    del qc_metrics

    # Evaluate all QC flags in one Polars operation
    df = evaluate_qc_flags(df, qc_flags)

    df, flags_df, output_columns = write_tabular_reports(df, qc_flags, accepted_tsv_columns, args.output_prefix)
    datavzrd_config = build_datavzrd_config(df, flags_df, output_columns, args.output_prefix)
    with open(f"{args.output_prefix}.qc_align_datavzrd.yml", 'w') as f:
        yaml.safe_dump(datavzrd_config, f, sort_keys=False)
    log(f"datavzrd config saved to {args.output_prefix}.qc_align_datavzrd.yml")


def write_tabular_reports(df, qc_flags, accepted_tsv_columns, output_prefix):
    """Write the two-header metric table, flag matrix and rule descriptions."""
    # Save the dict to a TSV file, retaining only columns defined in accepted_tsv_columns
    accepted_tsv_columns = accepted_tsv_columns + ['QC_warn_flags', 'QC_fail_flags', 'N_FAIL_conditions', 'N_WARN_conditions', 'QC_status']
    accepted = set(accepted_tsv_columns)
    df = df.select([col for col in df.columns if col in accepted])
    # Place QC_status right after the sample identifiers
    id_columns = [col for col in ('sample_id', 'bam_sample_id') if col in df.columns]
    df = df.select(id_columns + ['QC_status'] + [col for col in df.columns if col not in id_columns + ['QC_status']])

    # Polars does not support pandas-style MultiIndex headers, so write the two-row TSV header manually.
    output_columns = []
    for col in df.columns:
        if col in ('sample_id', 'bam_sample_id'):
            output_columns.append((col, 'sampleid_section'))
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
        'sample_id', 'bam_sample_id',
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
                'description': description
            })
    qc_criteria_df = pl.DataFrame(qc_criteria_rows, schema={
        'severity': pl.String, 'flag': pl.String, 'description': pl.String})
    qc_criteria_df.write_csv(f"{output_prefix}.qc_criteria.tsv", separator='\t')
    log(f"QC criteria file saved to {output_prefix}.qc_criteria.tsv")

    return df, qc_flags_status_df, output_columns


def build_datavzrd_config(df, qc_flags_status_df, output_columns, output_prefix):
    """Describe the existing report tables for datavzrd."""
    secondary_header = [section for _, section in output_columns]
    # Export a datavzrd config for the QC metric table. Columns are keyed by the primary header
    # (first row) and the secondary header (QC section) is rendered as a colored heatmap
    sections = list(dict.fromkeys(secondary_header))
    palette = ["#eeeeee", "#c6dbef", "#c7e9c0", "#fdd0a2", "#dadaeb", "#fcbba1", "#d9d9d9", "#bcbddc", "#fde0dd"]
    section_colors = [palette[i % len(palette)] for i in range(len(sections))]

    datavzrd_columns = {}
    for col, (primary, _) in zip(df.columns, output_columns):
        if col == 'sample_id':
            datavzrd_columns[primary] = {'label': 'Sample ID', 'display-mode': 'pinned'}
        elif col == 'bam_sample_id':
            datavzrd_columns[primary] = {'label': 'BAM Sample ID', 'display-mode': 'pinned'}
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
                    }
                }
        else:
            datavzrd_columns.setdefault(primary, {'label': primary, 'precision': 3})

    datavzrd_config = {
        'name': 'QC metrics report',
        'datasets': {
            'qc_metrics': {
                'path': f"{output_prefix}.tsv",
                'separator': '\t',
                'headers': 2,
            }
        },
        'views': {
            'qc_metrics': {
                'dataset': 'qc_metrics',
                'desc': 'Alignment QC metrics',
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
        'bam_sample_id': {'label': 'BAM Sample ID', 'display-mode': 'pinned'},
    }
    qc_flags_columns['QC_status'] = {
        'display-mode': 'pinned',
        'plot': {
            'heatmap': {
                'scale': 'ordinal',
                'domain': ['PASS', 'WARN', 'FAIL'],
                'range': ['#00a000', '#ffa500', '#ff0000'],
            }
        }
    }
    qc_flags_columns['N_WARN_conditions'] = {'display-mode': 'pinned'}
    qc_flags_columns['N_FAIL_conditions'] = {'display-mode': 'pinned'}
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
    datavzrd_config['datasets']['qc_flags'] = {
        'path': f"{output_prefix}.qc_flags.tsv",
        'separator': '\t',
    }
    datavzrd_config['views']['qc_flags'] = {
        'dataset': 'qc_flags',
        'desc': 'Alignment QC flags status (1 = flag raised)',
        'render-table': {'columns': qc_flags_columns},
    }

    # Add the QC criteria table: severity colored (fail_rules red, warning_rules orange)
    datavzrd_config['datasets']['qc_criteria'] = {
        'path': f"{output_prefix}.qc_criteria.tsv",
        'separator': '\t',
    }
    datavzrd_config['views']['qc_criteria'] = {
        'dataset': 'qc_criteria',
        'desc': 'Alignment QC criteria',
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
                    }
                },
                'flag': {'label': 'Flag'},
                'description': {'label': 'Description'},
            }
        },
    }

    return datavzrd_config


if __name__ == '__main__':
    main()
