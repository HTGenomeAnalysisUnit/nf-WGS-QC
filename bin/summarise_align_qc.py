"""Combine alignment QC manifests into a selected TSV and a detailed JSON report.

JSON stores nested reports per sample. 
TSV metrics use each parser's keys, joining nested keys with dash
"""

import os
import ast
import argparse
import io
import json
import operator
import tokenize
import yaml
import math
import polars as pl
import numpy as np
from datetime import datetime
from functools import reduce

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
            'reads MQ0',
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
            # 'dist_global-pct_bases_min_depth_{depth}' will be added downstream according to arguments
            # 'dist_region-pct_bases_min_depth_{depth}' will be added downstream according to arguments
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
        'visible': ['FREEMIX'],
        'hidden': []
    },
    'qrs' : {
        'visible': ['median_baseq', 'mean_read_gc_perc'],
        'hidden': ['n_runs', 'n_flowcells', 'n_lanes']
    }
}

def log(msg: str = "", level: str = "INFO"):
    print(f"[{datetime.now().isoformat()}] - {level} - {msg}")

def sanitize_key(key: str):
    # Change spaces into underscore, remove brackets and symbols
    clean_key = key.rstrip(':').translate(str.maketrans('', '', '()\\/')).replace(' ','_').lower()
    return clean_key.replace('%', 'pct').replace('&', 'and')

def csv_list(value: str) -> list[str]:
    items = [x.strip() for x in value.split(",") if x.strip()]
    if not items:
        raise argparse.ArgumentTypeError("Expected comma-separated values")
    return items

# A function to parse SN block from samtools stat output into a dictionary
# Parse samtools stats file into a dictionary with key-value pairs for each statistic
def parse_samtools_stats(file_path, base_q_thresholds=[30]):
    sn_block_done = False
    baseq_block_done = False
    stats = {}
    total_bases = { 'FFQ': 0, 'LFQ': 0}
    passing_q_bases = { 'FFQ': {}, 'LFQ': {} }
    for t in base_q_thresholds:
        passing_q_bases['FFQ'][t] = 0
        passing_q_bases['LFQ'][t] = 0
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            # Parse summary stats (SN) block
            if parts[0] == 'SN':
                key = sanitize_key(parts[1])
                value = float(parts[2].strip())
                stats[key] = value
                sn_block_done = True
            # Parse base quality blocks (FFQ and LFQ)
            elif parts[0] in ['FFQ', 'LFQ']:
                # parts[1] is the cycle; parts[2:] contains counts for Q0, Q1, ...
                section = parts[0]
                counts = list(map(int, parts[2:]))
                total_bases[section] += sum(counts)
                for t in base_q_thresholds:
                    passing_q_bases[section][t] += sum(counts[t:])
                baseq_block_done = True
            # Comments separate FFQ and LFQ; stop at the next data block.
            elif parts[0] and not line.startswith('#') and sn_block_done and baseq_block_done:
                break
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
    
    # Return the dict with all values rounded to precision set my PRECISION
    for key in stats:
        if isinstance(stats[key], float):
            stats[key] = round(stats[key], PRECISION)
    return stats

# Function to parse the samtools idxstats output into a dictionary
def parse_samtools_idxstats(file_path, autosomes, sex_chroms, canonical_chroms):
    autosomes_norm = [c.lstrip('chr') for c in autosomes]
    canonical_chroms_norm = [c.lstrip('chr') for c in canonical_chroms]
    idxstats = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                chrom_size = int(parts[1]) / 1e6
                mapped = int(parts[2])
                unmapped = int(parts[3])
                if (chrom.lstrip('chr') in canonical_chroms_norm + ['M', 'MT']):
                    idxstats[chrom] = {
                        'mapped': mapped, 
                        'mapped_norm': round(mapped / chrom_size if chrom_size > 0 else 0, PRECISION),
                        'unmapped': unmapped,
                        'fraction_mapped': round(mapped / (mapped + unmapped) if (mapped + unmapped) > 0 else 0, PRECISION)
                    }
    # Add collective stats for autosomes and canonical chromosomes collectively
    # Normalise idxstats keys and chromosome lists by stripping chr prefix for consistent lookup
    idxstats_norm = {k.lstrip('chr'): v for k, v in idxstats.items()}
    idxstats['autosomes'] = {'mapped': sum(idxstats_norm[chrom]['mapped'] for chrom in autosomes_norm if chrom in idxstats_norm),
                             'unmapped': sum(idxstats_norm[chrom]['unmapped'] for chrom in autosomes_norm if chrom in idxstats_norm),
                             'min_mapped_norm': min(idxstats_norm[chrom]['mapped_norm'] for chrom in autosomes_norm if chrom in idxstats_norm),
                             'min_mapped_norm_chrom': min((chrom for chrom in autosomes_norm if chrom in idxstats_norm), key=lambda chrom: idxstats_norm[chrom]['mapped_norm'] )}
    # Compute fraction of the min mapped_norm value relative to the mean for autosomes
    mean_value = np.mean([idxstats_norm[chrom]['mapped_norm'] for chrom in autosomes_norm if chrom in idxstats_norm])
    min_value = min([idxstats_norm[chrom]['mapped_norm'] for chrom in autosomes_norm if chrom in idxstats_norm])
    idxstats['autosomes']['lowest_cov_fraction_of_mean'] = round(float(min_value / mean_value if mean_value > 0 else 0), PRECISION)
    # Compute fraction mapped for autosomes
    total_reads_autosomes = idxstats['autosomes']['mapped'] + idxstats['autosomes']['unmapped']
    idxstats['autosomes']['fraction_mapped'] = round(idxstats['autosomes']['mapped'] / total_reads_autosomes if total_reads_autosomes > 0 else 0, PRECISION)

    # Repeat for all canonical chroms
    idxstats['canonical'] = {'mapped': sum(idxstats_norm[chrom]['mapped'] for chrom in canonical_chroms_norm if chrom in idxstats_norm),
                             'unmapped': sum(idxstats_norm[chrom]['unmapped'] for chrom in canonical_chroms_norm if chrom in idxstats_norm)}
    total_reads_canonical = idxstats['canonical']['mapped'] + idxstats['canonical']['unmapped']
    idxstats['canonical']['fraction_mapped'] = round(idxstats['canonical']['mapped'] / total_reads_canonical if total_reads_canonical > 0 else 0, PRECISION)
    
    # Now compite X and Y normalized coverage as mapped_norm / mapped_norm for an autosome with size similar to X
    sex_chrom_norm = [c.lstrip('chr') for c in sex_chroms]
    scaling_chrom = autosomes_norm[0]
    for c in sex_chrom_norm:
        chrom_label = c
        if f"chr{c}" in idxstats:
            chrom_label = f"chr{c}"
        idxstats[chrom_label]['scaled_cov'] = round((idxstats_norm[c]['mapped_norm'] / idxstats_norm[scaling_chrom]['mapped_norm']) * 2, PRECISION) if idxstats_norm[scaling_chrom]['mapped_norm'] > 0 else None
    
    return idxstats

# Function to parse generic reports made by header and a single line, like verifybamid selfSM
def parse_simple_report(file_path, columns=[]):
    header = []
    with open(file_path, 'r') as f:
        parsed_header = False
        for line in f:
            if not parsed_header:
                header = line.strip().split('\t')
                parsed_header = True
                continue
            values = dict(zip(header, line.strip().split('\t')))
        
    if len(columns) > 0:
        values = {k: v for k, v in values.items() if k in columns}
    
    # Sanitize keys
    result = {sanitize_key(k): v for k, v in values.items()}
    return result

# Function to parse mosdepth summary output into a dictionary
def parse_mosdepth_summary(file_path, autosomes, canonical_chroms):
    coverage = {}
    with open(file_path, 'r') as f:
        header = False
        for line in f:
            if not header:
                header = True
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            bases = int(parts[2])
            mean = float(parts[3])
            min_cov = int(parts[4])
            max_cov = int(parts[5])
            coverage[chrom] = {'bases': bases, 'mean': mean, 'min': min_cov, 'max': max_cov}
    # Similar as before now compute for autosomes and canonical chromosomes only
    coverage_norm = {k.lstrip('chr'): v for k, v in coverage.items()}
    autosomes_norm = [c.lstrip('chr') for c in autosomes]
    canonical_chroms_norm = [c.lstrip('chr') for c in canonical_chroms]
    coverage['autosomes'] = {'bases': sum(coverage_norm[chrom]['bases'] for chrom in autosomes_norm if chrom in coverage_norm),
                             'mean_cov': round(float(np.mean([coverage_norm[chrom]['mean'] for chrom in autosomes_norm if chrom in coverage_norm])), PRECISION),
                             'min_cov': min(coverage_norm[chrom]['min'] for chrom in autosomes_norm if chrom in coverage_norm),
                             'max_cov': max(coverage_norm[chrom]['max'] for chrom in autosomes_norm if chrom in coverage_norm)}
    coverage['canonical'] = {'bases': sum(coverage_norm[chrom]['bases'] for chrom in canonical_chroms_norm if chrom in coverage_norm),
                             'mean_cov': round(float(np.mean([coverage_norm[chrom]['mean'] for chrom in canonical_chroms_norm if chrom in coverage_norm])), PRECISION),
                             'min_cov': min(coverage_norm[chrom]['min'] for chrom in canonical_chroms_norm if chrom in coverage_norm),
                             'max_cov': max(coverage_norm[chrom]['max'] for chrom in canonical_chroms_norm if chrom in coverage_norm)}
    
    # Now compute coverage for all contigs that are not in the canonical set
    non_canonical_chroms_norm = [chrom for chrom in coverage_norm
                                if chrom not in canonical_chroms_norm + ['M', 'MT','total'] + [f"{c}_region" for c in canonical_chroms_norm + ['M', 'MT', 'total']]]
    coverage['extra_contigs'] = {'bases': sum(coverage_norm[chrom]['bases'] for chrom in non_canonical_chroms_norm if chrom in coverage_norm),
                                 'mean_cov': round(float(np.mean([coverage_norm[chrom]['mean'] for chrom in non_canonical_chroms_norm])), PRECISION) if non_canonical_chroms_norm else None,
                                 'min_cov': min((coverage_norm[chrom]['min'] for chrom in non_canonical_chroms_norm), default=None),
                                 'max_cov': max((coverage_norm[chrom]['max'] for chrom in non_canonical_chroms_norm), default=None)}

    # If some regions are present, compute coverage for autosomes and canonical chromosomes regions
    autosomes_norm_regions = [f"{chrom}_region" for chrom in autosomes_norm]
    canonical_chroms_norm_regions = [f"{chrom}_region" for chrom in canonical_chroms_norm]
    
    if not any(chrom in coverage_norm for chrom in autosomes_norm_regions + canonical_chroms_norm_regions):
        return coverage
    
    coverage['autosomes_region'] = {'bases': sum(coverage_norm[chrom]['bases'] for chrom in autosomes_norm_regions if chrom in coverage_norm),
                                      'mean_cov': round(float(np.mean([coverage_norm[chrom]['mean'] for chrom in autosomes_norm_regions if chrom in coverage_norm])), PRECISION),
                                      'min_cov': min(coverage_norm[chrom]['min'] for chrom in autosomes_norm_regions if chrom in coverage_norm),
                                      'max_cov': max(coverage_norm[chrom]['max'] for chrom in autosomes_norm_regions if chrom in coverage_norm)}
    coverage['canonical_region'] = {'bases': sum(coverage_norm[chrom]['bases'] for chrom in canonical_chroms_norm_regions if chrom in coverage_norm),
                                     'mean_cov': round(float(np.mean([coverage_norm[chrom]['mean'] for chrom in canonical_chroms_norm_regions if chrom in coverage_norm])), PRECISION),
                                     'min_cov': min(coverage_norm[chrom]['min'] for chrom in canonical_chroms_norm_regions if chrom in coverage_norm),
                                     'max_cov': max(coverage_norm[chrom]['max'] for chrom in canonical_chroms_norm_regions if chrom in coverage_norm)}
    
    # Remove the single extra contigs from coverage dict
    for chrom in non_canonical_chroms_norm:
        coverage.pop(chrom, None)
        coverage.pop(f"chr{chrom}", None)

    return coverage

# Function to parse mosdepth coverage distribution to a dictionary
def parse_mosdepth_distribution(file_path, target_depths, autosomes, canonical_chroms):
    coverage_distribution = {}
    with open(file_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            depth = int(fields[1])
            fraction = float(fields[2])
            if depth in target_depths:
                if chrom not in coverage_distribution:
                    coverage_distribution[chrom] = {}
                coverage_distribution[chrom][f"pct_bases_min_depth_{depth}"] = fraction

    # We also return mean fraction value across autosomes and canonical chromosomes collectively
    coverage_norm = {k.lstrip('chr'): v for k, v in coverage_distribution.items()}
    autosomes_norm = [c.lstrip('chr') for c in autosomes]
    canonical_chroms_norm = [c.lstrip('chr') for c in canonical_chroms]
    for depth in target_depths:
        autosomes_fractions = [coverage_norm[chrom][f"pct_bases_min_depth_{depth}"] for chrom in autosomes_norm if chrom in coverage_norm and f"pct_bases_min_depth_{depth}" in coverage_norm[chrom]]
        canonical_fractions = [coverage_norm[chrom][f"pct_bases_min_depth_{depth}"] for chrom in canonical_chroms_norm if chrom in coverage_norm and f"pct_bases_min_depth_{depth}" in coverage_norm[chrom]]
        if autosomes_fractions:
            coverage_distribution.setdefault('autosomes', {})[f"pct_bases_min_depth_{depth}"] = round(sum(autosomes_fractions) / len(autosomes_fractions), PRECISION)
        if canonical_fractions:
            coverage_distribution.setdefault('canonical', {})[f"pct_bases_min_depth_{depth}"] = round(sum(canonical_fractions) / len(canonical_fractions), PRECISION)
    
    # Remove the single extra contigs from coverage dict
    non_canonical_chroms_norm = [chrom for chrom in coverage_norm
                            if chrom not in canonical_chroms_norm + ['M', 'MT','total'] + [f"{c}_region" for c in canonical_chroms_norm + ['M', 'MT', 'total']]]
    for chrom in non_canonical_chroms_norm:
        coverage_distribution.pop(chrom, None)
        coverage_distribution.pop(f"chr{chrom}", None)
    
    return coverage_distribution

# Function to parse somalier pairs output file
def parse_somalier_pairs(file_path, relatedness_duplicate_value, relatedness_max_discrepancy, id_map={}):
    pairs_df = pl.read_csv(file_path, separator='\t')
    # Samples can be defined in either sample_a or sample_b columns and we want to have a dictionary that for each sample report the list of distinct sample with relatedness above the threshold
    related_samples = {}
    for row in pairs_df.iter_rows(named=True):
        sample_a = id_map.get(row['#sample_a'], row['#sample_a'])
        sample_b = id_map.get(row['sample_b'], row['sample_b'])
        relatedness = float(row['relatedness'])
        expected_relatedness = float(row['expected_relatedness'])
        related_samples.setdefault(sample_a, {"duplicated_samples": [], "unexpected_related": []})
        related_samples.setdefault(sample_b, {"duplicated_samples": [], "unexpected_related": []})
        # When relatedness is above the duplicate threshold, we consider the samples as duplicated
        if relatedness >= relatedness_duplicate_value:
            related_samples[sample_a]["duplicated_samples"].append(sample_b)
            related_samples[sample_b]["duplicated_samples"].append(sample_a)
        # When expected_relatedness is -1 and and relatedness is >= 0.05, or when the difference between relatedness and expected_relatedness is >= relatedness_max_discrepancy, we consider it unexpected relatedness
        if (expected_relatedness == -1 and relatedness >= 0.05) or (expected_relatedness != -1 and abs(relatedness - expected_relatedness) >= relatedness_max_discrepancy):
            related_samples[sample_a]["unexpected_related"].append(sample_b)
            related_samples[sample_b]["unexpected_related"].append(sample_a)

        # Convert to pipe separated strings for duplicated_samples and unexpected_related
        for sample in related_samples:
            related_samples[sample]["duplicated_samples"] = "|".join(related_samples[sample]["duplicated_samples"])
            related_samples[sample]["unexpected_related"] = "|".join(related_samples[sample]["unexpected_related"])

    return related_samples

# Function to parse somalier samples output file
def parse_somalier_samples(file_path, id_map={}):
    samples_df = pl.read_csv(file_path, separator='\t')
    samples = {}
    for row in samples_df.iter_rows(named=True):
        sample_id = id_map.get(row['sample_id'], row['sample_id'])
        samples[sample_id] = {
            'provided_sex': row['original_pedigree_sex'],
            'X_depth_mean': row['X_depth_mean'],
            'Y_depth_mean': row['Y_depth_mean'],
            'scaled_X': round((row['X_depth_mean'] / row['depth_mean']) * 2, PRECISION),
            'scaled_Y': round((row['Y_depth_mean'] / row['depth_mean']) * 2, PRECISION),
            'X_n_vars': row['X_n'],
            'X_het': row['X_het'],
            'Y_n_vars': row['Y_n']
        }
        inferred_sex = "male" if (
            (0.8 <= samples[sample_id]['scaled_Y'] <= 1.5) or
            (samples[sample_id]['X_het'] < 5 and samples[sample_id]['scaled_X'] <= 1.2)
            ) else "female"
        inferred_karyotype = infer_karyotype(samples[sample_id]['scaled_X'], samples[sample_id]['scaled_Y'], row['X_het'], row['Y_n'])
        samples[sample_id]['inferred_sex'] = inferred_sex
        samples[sample_id]['inferred_karyotype'] = inferred_karyotype
        samples[sample_id]['concordant_sex'] = 1 if (samples[sample_id]['provided_sex'].lower() == inferred_sex or samples[sample_id]['provided_sex'] == 'unknown') else 0
    return samples

# Function to parse somalier ancestry output file
def parse_somalier_ancestry(file_path, id_map={}):
    ancestry_df = pl.read_csv(file_path, separator='\t')
    # Filter the ancestry dataframe to keep only relevant samples, those where given_ancestry was empty
    ancestry_df = ancestry_df.filter(pl.col('given_ancestry').is_null() | (pl.col('given_ancestry') == ''))
    ancestry = {}
    for row in ancestry_df.iter_rows(named=True):
        sample_id = id_map.get(row['#sample_id'], row['#sample_id'])
        ancestry[sample_id] = {
            'ancestry': row['predicted_ancestry'],
            'ancestry_prob': row[f"{row['predicted_ancestry']}_prob"]
        }
    return ancestry

def infer_karyotype(scaled_X, scaled_Y, X_het, Y_n):
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

def read_table(file_path):
    # Keep identifiers such as 001 as strings in manifests and reports.
    return pl.read_csv(file_path, separator='\t', infer_schema=False)

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

def read_pedigree(file_path):
    pedigree = {}
    columns = ['family_id', 'individual_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
    with open(file_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 6:
                raise ValueError(f"Expected six PED columns in {file_path}: {line.strip()}")
            if fields[0].upper() == 'FID' and fields[1].upper() == 'IID':
                continue
            record = dict(zip(columns, fields[:6]))
            record['sex'] = int(record['sex'])
            sample_id = record['individual_id']
            if sample_id in pedigree and pedigree[sample_id] != record:
                raise ValueError(f"Conflicting PED entries for {sample_id}")
            pedigree[sample_id] = record
    return pedigree

def flatten_dict(values, prefix=''):
    flattened = {}
    for key, value in values.items():
        column = f'{prefix}-{key}' if prefix else str(key)
        if isinstance(value, dict):
            flattened.update(flatten_dict(value, column))
        else:
            flattened[column] = value
    return flattened

def json_values(value):
    # Missing numeric measurements must be JSON null rather than NaN/Infinity.
    if isinstance(value, dict):
        return {key: json_values(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_values(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value

def qc_condition_expr(condition, columns):
    """Translate Python conditions to Polars, accepting flattened column names.

    Known column names (including hyphenated names) take precedence over
    arithmetic. Strings are left intact. Calls and attribute access are rejected.
    """
    aliases = {}
    column_tokens = []
    for index, column in enumerate(columns):
        alias = f'__qc_column_{index}'
        aliases[alias] = column
        parts = [token.string for token in tokenize.generate_tokens(io.StringIO(column).readline)
                 if token.type not in (tokenize.NEWLINE, tokenize.ENDMARKER)]
        column_tokens.append((parts, alias))
    column_tokens.sort(key=lambda item: len(item[0]), reverse=True)

    tokens = list(tokenize.generate_tokens(io.StringIO(condition).readline))
    rewritten = []
    index = 0
    while index < len(tokens):
        for parts, alias in column_tokens:
            if [token.string for token in tokens[index:index + len(parts)]] == parts:
                rewritten.append((tokenize.NAME, alias))
                index += len(parts)
                break
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
    for severity in ('fail_rules', 'warning_rules'):
        for flag, condition in (qc_flags.get(severity) or {}).items():
            name = f'__qc_rule_{len(rules)}'
            try:
                expression = qc_condition_expr(condition, df.columns)
                result = df.with_columns(expression.alias(name)).get_column(name)
                if result.dtype != pl.Boolean:
                    raise ValueError('QC conditions must return booleans')
            except (ValueError, TypeError, SyntaxError, tokenize.TokenError, pl.exceptions.PolarsError) as error:
                raise ValueError(f'Invalid {severity} flag {flag!r} ({condition!r}): {error}') from error
            conditions.append(result.fill_null(False))
            rules.append((severity, flag, name))

    if not rules:
        return df.with_columns(
            pl.lit('').alias('QC_fail_flags'),
            pl.lit('').alias('QC_warn_flags'),
            pl.lit(0, dtype=pl.UInt32).alias('N_FAIL_conditions'),
            pl.lit(0, dtype=pl.UInt32).alias('N_WARN_conditions'),
            pl.lit('PASS').alias('QC_status'),
        )

    flags = pl.DataFrame(conditions).select(
        pl.concat_str([
            pl.when(pl.col(name)).then(pl.lit(flag)).otherwise(None)
            for group, flag, name in rules if group == 'fail_rules'
        ], separator='|', ignore_nulls=True).fill_null('').alias('QC_fail_flags'),
        pl.concat_str([
            pl.when(pl.col(name)).then(pl.lit(flag)).otherwise(None)
            for group, flag, name in rules if group == 'warning_rules'
        ], separator='|', ignore_nulls=True).fill_null('').alias('QC_warn_flags'),
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

def main(argv=None):
    # Arguments
    parser = argparse.ArgumentParser(description="Summarise alignment QC")
    parser.add_argument("--sample_id_col", help="Column name for sample ID in input files", default="id")
    parser.add_argument("--bam_sample_id_col", help="Column name for BAM sample ID in input files", default="bam_sample_id")
    parser.add_argument("--qc_flags", required=True, help="QC flags YAML file")
    #parser.add_argument("--input_ped", required=True, help="Input PED file")
    parser.add_argument("--somalier_ancestry", required=True, help="Somalier ancestry TSV report")
    parser.add_argument("--somalier_pairs", required=True, help="Somalier pairs TSV report")
    parser.add_argument("--relatedness_duplicate_value", type=float, default=0.85, help="Pairs with relatedness equal or greater to this value are considered duplicates")
    parser.add_argument("--relatedness_max_discrepancy", type=float, default=0.2, help="Maximum allowed discrepancy between computed and expected relatedness")
    parser.add_argument("--somalier_samples", required=True, help="Somalier samples TSV report")
    parser.add_argument("--input_mosdepth", required=True, help="Input mosdepth manifest file")
    parser.add_argument("--input_samtools_idxstats", required=True, help="Input samtools idxstats manifest file")
    parser.add_argument("--input_samtools_stats", required=True, help="Input samtools stats manifest file")
    parser.add_argument("--input_verifybamid", required=True, help="Input verifybamid manifest file")
    parser.add_argument("--input_qrs_sample_runs", required=True, help="Input QRS sample runs manifest file")
    parser.add_argument("--autosomes", type=csv_list, default=AUTOSOMES, help="Comma-separated list of autosomes")
    parser.add_argument("--sex_chroms", type=csv_list, default=SEX_CHROMS, help="Comma-separated list of sex chromosomes")
    parser.add_argument("--target_depths", type=csv_list, default=['10', '20', '30'], help="Comma-separated list of target depths to report (default: 10,20,30)")
    parser.add_argument("--target_base_qualities", type=csv_list, default=['30'], help="Comma-separated list of target bases qualities to report (default: 30)")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files (TSV, JSON, summary tables)")
    args = parser.parse_args(argv)
    try:
        target_depths = [int(value) for value in args.target_depths]
        target_base_qualities = [int(value) for value in args.target_base_qualities]
    except ValueError:
        parser.error('Target depths and base qualities must be comma-separated integers')
    if any(value < 0 for value in target_depths + target_base_qualities):
        parser.error('Target depths and base qualities must be non-negative')

    # Assemble list of canonical chromosomes
    autosomes = args.autosomes
    sex_chroms = args.sex_chroms
    canonical_chroms = autosomes + sex_chroms
    sample_id_col = args.sample_id_col
    bam_sample_id_col = args.bam_sample_id_col
    relatedness_duplicate_value = args.relatedness_duplicate_value
    relatedness_max_discrepancy = args.relatedness_max_discrepancy

    log(f"Sample ID column: {sample_id_col}")
    log(f"BAM Sample ID column: {bam_sample_id_col}")
    log(f"Canonical chromosomes: {canonical_chroms}")
    log(f"Relatedness duplicate value: {relatedness_duplicate_value}")
    log(f"Relatedness max discrepancy: {relatedness_max_discrepancy}")

    # Compile full list of keys accepted for the TSV report
    # These becomes column names as tool-column_name
    accepted_tsv_keys = TSV_SELECTED_KEYS.copy()
    for q in target_base_qualities:
        accepted_tsv_keys['samtools_stats']['visible'].append(f"bases_at_least_Q{q}")
        accepted_tsv_keys['samtools_stats']['hidden'].append(f"bases_at_least_Q{q}_FFQ")
        accepted_tsv_keys['samtools_stats']['hidden'].append(f"bases_at_least_Q{q}_LFQ")
    
    for d in target_depths:
        accepted_tsv_keys['mosdepth']['visible'].append(f"dist_global-pct_bases_min_depth_{d}")

    accepted_tsv_columns = ['sample_id', 'bam_sample_id']
    for tool, categories in accepted_tsv_keys.items():
        accepted_tsv_columns.extend([f"{tool}-{key}" for key in categories['visible'] + categories['hidden']])

    # Read QC flags from YAML file and log them
    with open(args.qc_flags, 'r') as f:
        qc_flags = yaml.safe_load(f) or {}
    log(f"QC flags: {qc_flags}")

    # Read the input files
    log("Reading input manifest files")
    manifests = {
        'samtools_stats': read_table(args.input_samtools_stats),
        'samtools_idxstats': read_table(args.input_samtools_idxstats),
        'mosdepth': read_table(args.input_mosdepth),
        'verifybamid': read_table(args.input_verifybamid),
        'qrs': read_table(args.input_qrs_sample_runs),
    }

    log("Finished reading input files")

    qc_metrics = {}
    id_map = {}
    for report, manifest in manifests.items():
        log(f"== {report.upper()} ==")
        for row in manifest.iter_rows(named=True):
            sample_id, bam_sample_id = row[sample_id_col], row[bam_sample_id_col]
            result = qc_metrics.get(sample_id, {'bam_sample_id': bam_sample_id})
            id_map[bam_sample_id] = sample_id
            input_file = row['output_file']
            log(f"Processing {report} for sample {sample_id}, input file {input_file}")
            if report == 'samtools_stats':
                result[report] = parse_samtools_stats(input_file, target_base_qualities)
            elif report == 'samtools_idxstats':
                result[report] = parse_samtools_idxstats(input_file, autosomes, sex_chroms, canonical_chroms)
            elif report == 'mosdepth':
                if input_file.endswith('.mosdepth.summary.txt'):
                    result.setdefault(report, {}).setdefault('summary',{}).update(parse_mosdepth_summary(input_file, autosomes, canonical_chroms))
                elif input_file.endswith('.mosdepth.global.dist.txt'):
                    result.setdefault(report, {}).setdefault('dist_global',{}).update(parse_mosdepth_distribution(input_file, target_depths, autosomes, canonical_chroms))
                elif input_file.endswith('.mosdepth.region.dist.txt'):
                    result.setdefault(report, {}).setdefault('dist_region',{}).update(parse_mosdepth_distribution(input_file, target_depths, autosomes, canonical_chroms))
                else:
                    raise ValueError(f"Unrecognised mosdepth report: {input_file}")
            elif report == 'verifybamid':
                result[report] = parse_simple_report(input_file, columns=['FREEMIX'])
            elif report == 'qrs':
                if input_file.endswith('.per_sample_run.tsv'):
                    result.setdefault(report, {}).update(parse_simple_report(input_file, columns=['N_RUNS', 'N_FLOWCELLS', 'N_LANES']))
                elif input_file.endswith('.read_stats.tsv'):
                    result.setdefault(report, {}).update(parse_simple_report(input_file, columns=['MEAN_BASEQ', 'MEDIAN_BASEQ', 'MEAN_READ_GC_PERC']))
                else:
                    log(f"Unrecognised QRS report: {input_file}", level="WARN")
            qc_metrics[sample_id] = result

    log("Read Somalier files")
    somalier_ancestry = parse_somalier_ancestry(args.somalier_ancestry,id_map)
    somalier_pairs = parse_somalier_pairs(args.somalier_pairs, relatedness_duplicate_value, relatedness_max_discrepancy, id_map)
    somalier_samples = parse_somalier_samples(args.somalier_samples, id_map)

    log("Finished reading all data")

    # Add the somalier dictionary keys to the corresponding samples in the results given that the sample ids should match
    for sample_id, result in somalier_ancestry.items():
        qc_metrics.setdefault(sample_id, {}).setdefault('somalier', {}).update(result)

    for sample_id, result in somalier_samples.items():
        qc_metrics.setdefault(sample_id, {}).setdefault('somalier', {}).update(result)

    for sample_id, result in somalier_pairs.items():
        qc_metrics.setdefault(sample_id, {}).setdefault('somalier', {}).update(result)

    # Save the full stats in qc_metrics[sample_id]: dict() to JSON 
    with open(f"{args.output_prefix}.json", 'w') as f:
        json.dump(qc_metrics, f, indent=4)
    log(f"JSON file saved to {args.output_prefix}.json")

    # Flatten the nested dictionaries 
    flattened_qc_metrics = [
        {'sample_id': sample_id, **flatten_dict(metrics)}
        for sample_id, metrics in qc_metrics.items()
    ]

    # Convert the flattened QC metrics to a Polars DataFrame
    df = pl.DataFrame(flattened_qc_metrics)

    # Evaluate QC f;ags
    df = evaluate_qc_flags(df, qc_flags)
    
    # Save the dict to a TSV file, retaining only columns defined in accepted_tsv_columns
    accepted_tsv_columns += ['QC_warn_flags', 'QC_fail_flags', 'N_FAIL_conditions', 'N_WARN_conditions', 'QC_status']
    df = df.select([col for col in df.columns if col in accepted_tsv_columns])
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
    output_path = f"{args.output_prefix}.tsv"
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
    log(f"TSV file saved to {args.output_prefix}.tsv")

    # Save a table reporting for each sample_id, bam_sample_id the 0/1 status for each FAIL and WARN QC flag, the number of warn / fail active flagas and the final QC status
    flag_names_warn = dict.fromkeys(
        flag for severity in ('warning_rules',)
        for flag in (qc_flags.get(severity) or {})
    )
    flag_names_fail = dict.fromkeys(
        flag for severity in ('fail_rules',)
        for flag in (qc_flags.get(severity) or {})
    )
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
    qc_flags_status_df.write_csv(f"{args.output_prefix}.qc_flags.tsv", separator='\t')
    log(f"QC flags status file saved to {args.output_prefix}.qc_flags.tsv")

    # Dump the QC criteria YML as a TSV table with columns for severity, flag name and description
    qc_criteria_rows = []
    for severity, flags in qc_flags.items():
        for flag, description in flags.items():
            qc_criteria_rows.append({
                'severity': severity,
                'flag': flag,
                'description': description
            })
    qc_criteria_df = pl.DataFrame(qc_criteria_rows)
    qc_criteria_df.write_csv(f"{args.output_prefix}.qc_criteria.tsv", separator='\t')
    log(f"QC criteria file saved to {args.output_prefix}.qc_criteria.tsv")

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
                'path': f"{args.output_prefix}.tsv",
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
        'path': f"{args.output_prefix}.qc_flags.tsv",
        'separator': '\t',
    }
    datavzrd_config['views']['qc_flags'] = {
        'dataset': 'qc_flags',
        'desc': 'Alignment QC flags status (1 = flag raised)',
        'render-table': {'columns': qc_flags_columns},
    }

    # Add the QC criteria table: severity colored (fail_rules red, warning_rules orange)
    datavzrd_config['datasets']['qc_criteria'] = {
        'path': f"{args.output_prefix}.qc_criteria.tsv",
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

    with open(f"{args.output_prefix}.qc_align_datavzrd.yml", 'w') as f:
        yaml.safe_dump(datavzrd_config, f, sort_keys=False)
    log(f"datavzrd config saved to {args.output_prefix}.qc_align_datavzrd.yml")

if __name__ == '__main__':
    main()
