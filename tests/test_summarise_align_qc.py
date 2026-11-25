"""Regression tests using generated inputs; run with unittest discovery."""

import contextlib
import importlib.util
import io
import itertools
import json
from pathlib import Path
import tempfile
import unittest

import polars as pl
import yaml

SCRIPT = Path(__file__).resolve().parents[1] / 'bin' / 'summarise_align_qc.py'
SPEC = importlib.util.spec_from_file_location('summarise_align_qc', SCRIPT)
qc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(qc)


def make_cohort(root, samples=3, cycles=2, max_depth=40):
    """Generate eight metric files per sample plus cohort-level Somalier files."""
    manifests = {name: ['id\tbam_sample_id\toutput_file'] for name in
                 ('samtools_stats', 'samtools_idxstats', 'mosdepth', 'verifybamid', 'qrs_sample_runs')}
    ids = [f'{index:05d}' for index in range(samples)]
    somalier_samples = ['sample_id\toriginal_pedigree_sex\tX_depth_mean\tY_depth_mean\tdepth_mean\tX_n\tX_het\tY_n']
    ancestry = ['#sample_id\tgiven_ancestry\tpredicted_ancestry\tEUR_prob']
    stats = ''.join(f'SN\t{key}:\t{value}\n' for key, value in {
        'sequences': 1000, 'reads duplicated': 100, 'reads mapped': 900,
        'reads properly paired': 800, 'average quality': 36, 'average length': 150,
        'total length': 150000, 'error rate': 0.001,
    }.items())
    for section in ('FFQ', 'LFQ'):
        stats += '# quality counts\n'
        stats += ''.join(f'{section}\t{cycle}\t' + '\t'.join(str(q) for q in range(41)) + '\n'
                         for cycle in range(1, cycles + 1))
    stats += 'IS\t0\t0\n'
    idxstats = ''.join(f'chr{c}\t1000000\t100\t0\n' for c in ('1', '2', 'X', 'Y'))
    idxstats += '*\t0\t0\t50\n'
    summary = 'chrom\tlength\tbases\tmean\tmin\tmax\n'
    summary += ''.join(f'chr{c}{suffix}\t100\t2500\t25\t0\t50\n'
                       for c in ('1', '2', 'X', 'Y') for suffix in ('', '_region'))
    summary += 'total\t400\t10000\t25\t0\t50\ntotal_region\t400\t10000\t25\t0\t50\n'
    summary += 'chrUn_gl\t10\t20\t2\t0\t5\n'
    distribution = ''.join(f'{chrom}\t{depth}\t{max(0, 1-depth/100):.4f}\n'
                           for chrom in ('chr1', 'chr2', 'chrX', 'chrY', 'total', 'chrUn_gl')
                           for depth in range(max_depth, -1, -1))
    reports = {
        'stats': ('samtools_stats', stats),
        'idxstats': ('samtools_idxstats', idxstats),
        'mosdepth.summary.txt': ('mosdepth', summary),
        'mosdepth.global.dist.txt': ('mosdepth', distribution),
        'mosdepth.region.dist.txt': ('mosdepth', distribution),
        'selfSM': ('verifybamid', 'FREEMIX\tUNUSED\n0.005\tignored\n'),
        'per_sample_run.tsv': ('qrs_sample_runs', 'N_RUNS\tN_FLOWCELLS\tN_LANES\n1\t2\t3\n'),
        'read_stats.tsv': ('qrs_sample_runs', 'MEAN_BASEQ\tMEDIAN_BASEQ\tMEAN_READ_GC_PERC\n35\t37\t0.41\n'),
    }
    for sample in ids:
        for suffix, (manifest, contents) in reports.items():
            path = root / f'{sample}.{suffix}'
            path.write_text(contents)
            manifests[manifest].append(f'sample_{sample}\t{sample}\t{path}')
        somalier_samples.append(f'{sample}\tunknown\t10\t10\t20\t100\t0\t20')
        ancestry.append(f'{sample}\t\tEUR\t0.99')
    args = ['--autosomes', '1,2', '--sex_chroms', 'X,Y']
    for name, lines in manifests.items():
        path = root / f'{name}.tsv'
        path.write_text('\n'.join(lines) + '\n')
        args.extend([f'--input_{name}', str(path)])
    pairs = ['#sample_a\tsample_b\trelatedness\texpected_relatedness']
    pairs.extend(f'{a}\t{b}\t0\t-1' for a, b in itertools.combinations(ids, 2))
    for name, lines in [('samples', somalier_samples), ('ancestry', ancestry), ('pairs', pairs)]:
        path = root / f'somalier_{name}.tsv'
        path.write_text('\n'.join(lines) + '\n')
        args.extend([f'--somalier_{name}', str(path)])
    flags = root / 'flags.yml'
    flags.write_text(yaml.safe_dump({
        'fail_rules': {'low_coverage': 'mosdepth-summary-autosomes-mean_cov < 20'},
        'warning_rules': {'duplicate': "somalier-duplicated_samples != ''"},
    }))
    return args + ['--qc_flags', str(flags), '--output_prefix', str(root / 'report')]


class SummariseAlignQCTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write(self, name, contents):
        path = self.root / name
        path.write_text(contents)
        return path

    def run_main(self, args):
        with contextlib.redirect_stdout(io.StringIO()):
            qc.main(args)

    def test_pairs_multiple_flagged_rows_and_string_ids(self):
        path = self.write('pairs.tsv', '#sample_a\tsample_b\trelatedness\texpected_relatedness\n'
                          '001\t002\t0.9\t-1\n001\t003\t0.3\t0\n'
                          '001\t002\t0.9\t-1\n002\t003\t0\t-1\n'
                          '004\t005\t0.05\t-1\n006\t007\t0.25\t0.25\n')
        result = qc.parse_somalier_pairs(path, 0.85, 0.2, {'001': 'user_001'})
        self.assertEqual(result['user_001']['duplicated_samples'], '002')
        self.assertEqual(result['user_001']['unexpected_related'], '002|003')
        self.assertEqual(result['003']['unexpected_related'], 'user_001')
        self.assertEqual(result['004']['unexpected_related'], '005')
        self.assertEqual(result['006']['unexpected_related'], '')
        self.assertEqual(result['007']['duplicated_samples'], '')

    def test_pairs_header_only(self):
        path = self.write('pairs.tsv', '#sample_a\tsample_b\trelatedness\texpected_relatedness\n')
        self.assertEqual(qc.parse_somalier_pairs(path, 0.85, 0.2), {})

    def test_samtools_quality_thresholds_and_zero_counts(self):
        path = self.write('stats', 'SN\tsequences:\t10\nSN\treads duplicated:\t2\n'
                          'SN\treads mapped:\t8\nSN\treads properly paired:\t6\n'
                          'FFQ\t1\t1\t2\t3\n# LFQ follows\nLFQ\t1\t3\t2\t1\nIS\t0\t0\n')
        result = qc.parse_samtools_stats(path, (0, 1, 2, 30))
        self.assertEqual(result['bases_at_least_Q0'], 1)
        self.assertEqual(result['bases_at_least_Q1'], round(8/12, qc.PRECISION))
        self.assertEqual(result['bases_at_least_Q2_FFQ'], 0.5)
        self.assertEqual(result['bases_at_least_Q30'], 0)
        self.assertEqual(result['reads_pct_dup_on_mapped'], 0.25)
        zero = self.write('zero.stats', 'SN\tsequences:\t0\nSN\treads duplicated:\t0\n'
                          'SN\treads mapped:\t0\nSN\treads properly paired:\t0\n')
        self.assertIsNone(qc.parse_samtools_stats(zero)['reads_pct_mapped'])
        qc.write_metrics_json({'a': qc.parse_samtools_stats(zero)}, self.root/'zero.json')
        data = json.loads((self.root/'zero.json').read_text(), parse_constant=lambda s: self.fail(s))
        self.assertIsNone(data['a']['bases_at_least_Q30'])

    def test_distribution_filters_contigs_and_averages_available_values(self):
        path = self.write('dist', 'chr1\t10\t0.8\nchr2\t10\t0.6\nchrX\t10\t0.4\n'
                          'total\t10\t0.7\nchrUn\t10\t0.9\nchr1\t11\t0.75\nchr2\t20\t0.3\n')
        result = qc.parse_mosdepth_distribution(path, [10, 20], ['chr1', '2'], ['1', '2', 'X'])
        self.assertEqual(result['autosomes']['pct_bases_min_depth_10'], 0.7)
        self.assertEqual(result['canonical']['pct_bases_min_depth_10'], 0.6)
        self.assertEqual(result['autosomes']['pct_bases_min_depth_20'], 0.3)
        self.assertNotIn('chrUn', result)
        self.assertNotIn('pct_bases_min_depth_11', result['chr1'])

    def test_summary_keeps_unweighted_means_and_exact_contig_names(self):
        path = self.write('summary', 'chrom\tlength\tbases\tmean\tmin\tmax\n'
                          'chr1\t10\t100\t10\t0\t20\nchr2\t100\t3000\t30\t1\t40\n'
                          'chr1_region\t5\t50\t10\t0\t10\n'
                          'chr2_region\t5\t150\t30\t0\t30\nchrchrUn\t5\t25\t5\t0\t5\n')
        result = qc.parse_mosdepth_summary(path, ['1', '2'], ['1', '2'])
        self.assertEqual(result['autosomes']['mean_cov'], 20)
        self.assertEqual(result['extra_contigs']['bases'], 25)
        self.assertNotIn('chrchrUn', result)
        self.assertEqual(qc.normalize_chrom('chrchrUn'), 'chrUn')

    def test_idxstats_missing_sex_chroms_and_zero_autosomes(self):
        path = self.write('idxstats', 'chr1\t1000\t0\t0\n*\t0\t0\t10\n')
        result = qc.parse_samtools_idxstats(path, ['1'], ['X', 'Y'], ['1', 'X', 'Y'])
        self.assertEqual(result['autosomes']['min_mapped_norm'], 0)
        self.assertEqual(result['canonical']['fraction_mapped'], 0)

    def test_simple_report_numeric_values_and_empty_input(self):
        path = self.write('selfSM', 'FREEMIX\tUNUSED\n0.125\ttext\n')
        self.assertEqual(qc.parse_simple_report(path, ['FREEMIX']), {'freemix': 0.125})
        with self.assertRaisesRegex(ValueError, 'No data rows'):
            qc.parse_simple_report(self.write('empty', 'FREEMIX\n'))

    def test_zero_depth_somalier_sample(self):
        path = self.write('samples', 'sample_id\toriginal_pedigree_sex\tX_depth_mean\tY_depth_mean\tdepth_mean\tX_n\tX_het\tY_n\n'
                          '001\tunknown\t0\t0\t0\t0\t0\t0\n')
        result = qc.parse_somalier_samples(path)['001']
        self.assertIsNone(result['scaled_X'])
        self.assertEqual(result['inferred_karyotype'], 'unknown')

    def test_qc_longest_names_arithmetic_literals_and_nulls(self):
        df = pl.DataFrame({'a': [1, 2, None], 'a-b': [3, 1, None],
                           'value': ['XX', 'XY', 'NA']})
        result = qc.evaluate_qc_flags(df, {
            'fail_rules': {'numeric': 'a-b < 2 and a + 1 >= 3'},
            'warning_rules': {'missing': "value == 'NA'", 'sex': "value not in ['XX', 'XY', 'NA']"},
        })
        self.assertEqual(result['QC_status'].to_list(), ['PASS', 'FAIL', 'WARN'])
        self.assertEqual(result['N_FAIL_conditions'].to_list(), [0, 1, 0])
        self.assertEqual(result['QC_warn_flags'].to_list(), ['', '', 'missing'])

    def test_qc_one_severity_no_rules_and_scalar_conditions(self):
        df = pl.DataFrame({'a': [1, 2, 3]})
        for rules, status in [({}, 'PASS'), ({'fail_rules': {'always': 'True'}}, 'FAIL'),
                              ({'warning_rules': {'always': 'True'}}, 'WARN')]:
            with self.subTest(rules=rules):
                result = qc.evaluate_qc_flags(df, rules)
                self.assertEqual(result['QC_status'].to_list(), [status]*3)

    def test_qc_rejects_calls_and_unknown_columns(self):
        df = pl.DataFrame({'a': [1]})
        for condition in ('abs(a) > 1', '__import__("os")', 'missing > 1', 'a + 1'):
            with self.subTest(condition=condition), self.assertRaises(ValueError):
                qc.evaluate_qc_flags(df, {'fail_rules': {'bad': condition}})

    def test_frame_late_metrics_and_columns_used_only_by_rules(self):
        metrics = {str(i): {'bam_sample_id': str(i), 'metric': {'mean': 10}} for i in range(150)}
        metrics['149']['metric']['late'] = 0.5
        metrics['149']['unused'] = {'detail': 123}
        flags = {'fail_rules': {'late': 'metric-late > 0.1'}}
        df = qc.build_metrics_frame(metrics, ['metric-mean'], flags)
        self.assertNotIn('unused-detail', df.columns)
        result = qc.evaluate_qc_flags(df, flags)
        self.assertEqual(result['QC_status'].to_list(), ['PASS']*149 + ['FAIL'])

    def test_nonfinite_measurements_do_not_trigger_numeric_flags(self):
        metrics = {'001': {'value': float('nan')}, '002': {'value': float('inf')}}
        flags = {'fail_rules': {'high': 'value > 1'}}
        result = qc.evaluate_qc_flags(qc.build_metrics_frame(metrics, ['value'], flags), flags)
        self.assertEqual(result['QC_status'].to_list(), ['PASS', 'PASS'])

    def test_duplicate_manifest_rows_last_report_wins(self):
        args = make_cohort(self.root, samples=1)
        path = self.write('updated.selfSM', 'FREEMIX\n0.1\n')
        with (self.root/'verifybamid.tsv').open('a') as handle:
            handle.write(f'sample_00000\t00000\t{path}\n')
        self.run_main(args + ['--workers', '2'])
        data = json.loads((self.root/'report.json').read_text())
        self.assertEqual(data['sample_00000']['verifybamid']['freemix'], 0.1)

    def test_worker_parse_error_identifies_sample_and_file(self):
        args = make_cohort(self.root, samples=1)
        path = self.root/'00000.selfSM'
        path.write_text('FREEMIX\n')
        with self.assertRaisesRegex(ValueError, r"verifybamid.*sample_00000.*00000.selfSM.*No data rows"):
            self.run_main(args + ['--workers', '2'])

    def test_complete_reports_workers_repeated_runs_and_numeric_flags(self):
        args = make_cohort(self.root)
        flags = yaml.safe_load((self.root/'flags.yml').read_text())
        flags['fail_rules']['contamination'] = 'verifybamid-freemix > 0.01'
        (self.root/'flags.yml').write_text(yaml.safe_dump(flags))
        self.run_main(args)
        suffixes = ['json', 'tsv', 'qc_flags.tsv', 'qc_criteria.tsv', 'qc_align_datavzrd.yml']
        first = {suffix: (self.root/f'report.{suffix}').read_bytes() for suffix in suffixes}
        self.run_main(args + ['--workers', '3'])
        self.assertEqual(first, {suffix: (self.root/f'report.{suffix}').read_bytes() for suffix in suffixes})
        self.assertEqual(qc.TSV_SELECTED_KEYS['samtools_stats']['visible'].count('bases_at_least_Q30'), 0)
        data = json.loads(first['json'])
        self.assertEqual(data['sample_00000']['bam_sample_id'], '00000')
        self.assertEqual(data['sample_00000']['somalier']['ancestry'], 'EUR')
        self.assertEqual(data['sample_00000']['verifybamid']['freemix'], 0.005)
        header = first['tsv'].decode().splitlines()[0].split('\t')
        self.assertEqual(header[:3], ['sample_id', 'bam_sample_id', 'QC_status'])
        self.assertIn('freemix', header)
        self.assertIn('dist_region-total-pct_bases_min_depth_10', header)
        table = pl.read_csv(self.root/'report.qc_flags.tsv', separator='\t')
        self.assertEqual(table['QC_status'].to_list(), ['PASS']*3)
        config = yaml.safe_load(first['qc_align_datavzrd.yml'])
        self.assertEqual(config['datasets']['qc_metrics']['headers'], 2)

    def test_complete_report_without_qc_rules(self):
        args = make_cohort(self.root)
        (self.root/'flags.yml').write_text('{}\n')
        self.run_main(args)
        self.assertEqual((self.root/'report.qc_criteria.tsv').read_text(), 'severity\tflag\tdescription\n')


if __name__ == '__main__':
    unittest.main()
