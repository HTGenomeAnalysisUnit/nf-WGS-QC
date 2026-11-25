"""Regression tests for summarise_var_qc.py; run with unittest discovery."""

import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest

import polars as pl
import yaml

SCRIPT = Path(__file__).resolve().parents[1] / 'bin' / 'summarise_var_qc.py'
SPEC = importlib.util.spec_from_file_location('summarise_var_qc', SCRIPT)
qc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(qc)


def make_cohort(root, samples=3, include_regions=True):
    """Generate synthetic bcftools idxstats, bcftools stats (global + regions), and CHARR reports."""
    manifest_names = ['bcftools_idxstats', 'bcftools_stats', 'charr']
    if include_regions:
        manifest_names.append('bcftools_stats_regions')
    manifests = {name: ['id\tvcf_sample_id\toutput_file'] for name in manifest_names}

    ids = [f'{index:05d}' for index in range(samples)]
    idxstats = (
        'chr1\t2000000\t400\n'
        'chr2\t1000000\t150\n'
        'chrX\t1000000\t80\n'
        'chrY\t500000\t10\n'
        'chrM\t16569\t5\n'
        'chrUn_gl000220v1\t100000\t25\n'
    )
    stats = (
        '# bcftools stats\n'
        'ID\t0\tsample.bcf\n'
        'SN\t0\tnumber of samples:\t1\n'
        'SN\t0\tnumber of records:\t670\n'
        'SN\t0\tnumber of no-ALTs:\t0\n'
        'SN\t0\tnumber of SNPs:\t500\n'
        'SN\t0\tnumber of MNPs:\t0\n'
        'SN\t0\tnumber of indels:\t170\n'
        'SN\t0\tnumber of others:\t0\n'
        'SN\t0\tnumber of multiallelic sites:\t12\n'
        'SN\t0\tnumber of multiallelic SNP sites:\t3\n'
        'TSTV\t0\t340\t160\t2.13\t338\t160\t2.11\n'
        'QUAL\t0\t5.0\t50\t30\t20\t20\n'
        'QUAL\t0\t15.0\t100\t70\t30\t30\n'
        'QUAL\t0\t20.0\t50\t35\t15\t20\n'
        'QUAL\t0\t25.5\t300\t205\t95\t100\n'
        'ST\t0\tA>C\t20\n'
        'ST\t0\tA>G\t170\n'
        'ST\t0\tC>T\t170\n'
        'ST\t0\tG>T\t140\n'
        'DP\t0\t5\t70\t10.447761\t0\t0.0\n'
        'DP\t0\t10\t100\t14.925373\t0\t0.0\n'
        'DP\t0\t25\t400\t59.701493\t0\t0.0\n'
        'DP\t0\t120\t80\t11.940299\t0\t0.0\n'
        'DP\t0\t>500\t20\t2.985075\t0\t0.0\n'
    )
    region_stats = stats.replace('number of SNPs:\t500', 'number of SNPs:\t420')

    for i, sample in enumerate(ids):
        charr_val = 0.03 if i == samples - 1 else 0.005
        charr = (
            '## 500 variants considered for analysis\n'
            '#SAMPLE\tHOM_TOTAL\tHOM_COVERED\tHQ_HOM\tHQ_HOM_RATE\tHET_TOTAL\tHET_COVERED\tHQ_HET\tHQ_HET_RATE\tCHARR\tMEAN_REF_AB_HOM_ALT\tHETEROZYGOSITY_RATE\tINCONSISTENT_AB_HET_RATE\n'
            f'{sample}\t200\t190\t188\t0.98947\t300\t290\t285\t0.98276\t{charr_val}\t0.0015\t0.6\t0.012\n'
        )
        file_specs = {
            'bcftools_idxstats': ('bcftools_idxstats.txt', idxstats),
            'bcftools_stats': ('bcftools_stats.txt', stats),
            'charr': ('charr.tsv', charr),
        }
        if include_regions:
            file_specs['bcftools_stats_regions'] = ('regions.bcftools_stats.txt', region_stats)

        for manifest_key, (suffix, contents) in file_specs.items():
            path = root / f'{sample}.{suffix}'
            path.write_text(contents)
            manifests[manifest_key].append(f'user_{sample}\t{sample}\t{path}')

    args = ['--autosomes', '1,2', '--sex_chroms', 'X,Y']
    for name, lines in manifests.items():
        path = root / f'{name}.tsv'
        path.write_text('\n'.join(lines) + '\n')
        args.extend([f'--input_{name}', str(path)])

    flags = root / 'var_flags.yml'
    flags.write_text(yaml.safe_dump({
        'fail_rules': {'high_charr': 'charr-charr > 0.02'},
        'warning_rules': {'low_autosome_vars': 'bcftools_idxstats-autosomes-min_n_vars < 100'},
    }))
    return args + ['--qc_flags', str(flags), '--output_prefix', str(root / 'var_report')]


class SummariseVarQCTests(unittest.TestCase):
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

    def test_bcftools_idxstats_counts_normalization_and_missing_autosomes(self):
        path = self.write(
            'idxstats.txt',
            'chr1\t2000000\t400\n'
            'chr2\t1000000\t150\n'
            'chrX\t1000000\t80\n'
            'chrM\t16569\t5\n'
            'chrUn_1\t50000\t25\n',
        )
        result = qc.parse_bcftools_idxstats(path, ['1', '2', '3'], ['1', '2', '3', 'X', 'Y'])
        self.assertEqual(result['chr1'], {'n_vars': 400, 'n_vars_norm': 200.0})
        self.assertEqual(result['chr2'], {'n_vars': 150, 'n_vars_norm': 150.0})
        # chr3 and chrY were absent from bcftools index --stats -> filled with 0
        self.assertEqual(result['chr3'], {'n_vars': 0, 'n_vars_norm': 0})
        self.assertEqual(result['chrY'], {'n_vars': 0, 'n_vars_norm': 0})
        self.assertEqual(result['autosomes']['n_vars'], 550)
        self.assertEqual(result['autosomes']['min_n_vars'], 0)
        self.assertEqual(result['autosomes']['min_n_vars_chrom'], '3')
        self.assertEqual(result['autosomes']['mean_n_vars_norm'], round(350 / 3, qc.PRECISION))
        self.assertEqual(result['canonical']['n_vars'], 630)
        self.assertEqual(result['extra_contigs']['n_vars'], 25)
        self.assertEqual(result['total']['n_vars'], 660)

    def test_charr_comments_matching_sample_and_nan(self):
        path = self.write(
            'charr.tsv',
            '## comment line\n'
            '#SAMPLE\tCHARR\tHETEROZYGOSITY_RATE\tINCONSISTENT_AB_HET_RATE\tHQ_HOM_RATE\n'
            'S1\t0.005\t0.61\t0.01\tNaN\n'
            'S2\t0.025\t0.65\t0.04\t0.99\n',
        )
        s1 = qc.parse_charr_report(path, vcf_sample_id='S1')
        self.assertEqual(s1['charr'], 0.005)
        self.assertEqual(s1['heterozygosity_rate'], 0.61)
        self.assertEqual(s1['inconsistent_ab_het_rate'], 0.01)
        self.assertIsNone(s1['hq_hom_rate'])
        self.assertNotIn('sample', s1)

        s2 = qc.parse_charr_report(path, vcf_sample_id='S2')
        self.assertEqual(s2['charr'], 0.025)
        with self.assertRaisesRegex(ValueError, 'not found'):
            qc.parse_charr_report(path, vcf_sample_id='S3')

    def test_bcftools_stats_qual_and_dp_bins(self):
        path = self.write(
            'stats.txt',
            'SN\t0\tnumber of SNPs:\t500\n'
            'SN\t0\tnumber of indels:\t170\n'
            'SN\t0\tnumber of multiallelic sites:\t12\n'
            'TSTV\t0\t340\t160\t2.13\t338\t160\t2.11\n'
            'QUAL\t0\t.\t10\t5\t5\t0\n'
            'QUAL\t0\t10.0\t90\t60\t30\t20\n'
            'QUAL\t0\t20.0\t100\t65\t35\t50\n'
            'QUAL\t0\t35.2\t300\t210\t90\t100\n'
            'ST\t0\tA>G\t170\n'
            'ST\t0\tC>T\t170\n'
            'DP\t0\t0\t50\t5.0\t0\t0.0\n'
            'DP\t0\t10\t150\t15.0\t0\t0.0\n'
            'DP\t0\t50\t600\t60.0\t0\t0.0\n'
            'DP\t0\t100\t100\t10.0\t0\t0.0\n'
            'DP\t0\t150\t70\t7.0\t0\t0.0\n'
            'DP\t0\t>500\t30\t3.0\t0\t0.0\n',
        )
        res = qc.parse_bcftools_stats(path, var_qual_report=20, var_qual_dp=10)
        self.assertEqual(res['sn']['number_of_snps'], 500)
        self.assertEqual(res['tstv']['ts_tv'], round(340 / 160, qc.PRECISION))
        self.assertEqual(res['st'], {'A>G': 170, 'C>T': 170})

        # QUAL bins: 0, 10, 20, 30 (max_finite_qual = 35.2)
        # Missing QUAL ('.') has 10 SNPs that do not pass > 0
        self.assertEqual(res['qual']['n_snps_gt_0'], 490)
        self.assertEqual(res['qual']['fraction_snps_gt_0'], round(490 / 500, qc.PRECISION))
        self.assertEqual(res['qual']['n_snps_gt_10'], 400)
        self.assertEqual(res['qual']['n_snps_gt_20'], 300)
        self.assertEqual(res['qual']['fraction_snps_gt_20'], 0.6)
        self.assertEqual(res['qual']['n_indels_gt_20'], 100)
        self.assertEqual(res['qual']['fraction_indels_gt_20'], round(100 / 170, qc.PRECISION))
        self.assertEqual(res['qual']['n_snps_gt_30'], 300)
        self.assertNotIn('n_snps_gt_40', res['qual'])

        # DP bins: 0, 10, ..., 100 (total = 1000; >100 combines 150 and >500 = 100)
        self.assertEqual(res['dp']['n_genotypes_gt_0'], 950)
        self.assertEqual(res['dp']['fraction_genotypes_gt_0'], 0.95)
        self.assertEqual(res['dp']['n_genotypes_gt_10'], 800)
        self.assertEqual(res['dp']['fraction_genotypes_gt_10'], 0.8)
        self.assertEqual(res['dp']['n_genotypes_gt_90'], 200)
        self.assertEqual(res['dp']['n_genotypes_gt_100'], 100)
        self.assertEqual(res['dp']['fraction_genotypes_gt_100'], 0.1)

    def test_end_to_end_with_and_without_regions_and_workers(self):
        args = make_cohort(self.root, samples=3, include_regions=True)
        self.run_main(args)
        suffixes = ['json', 'tsv', 'qc_flags.tsv', 'qc_criteria.tsv', 'qc_var_datavzrd.yml']
        first = {suffix: (self.root / f'var_report.{suffix}').read_bytes() for suffix in suffixes}
        self.run_main(args + ['--workers', '3'])
        self.assertEqual(first, {suffix: (self.root / f'var_report.{suffix}').read_bytes() for suffix in suffixes})

        data = json.loads(first['json'])
        self.assertEqual(data['user_00000']['vcf_sample_id'], '00000')
        self.assertIn('bcftools_stats_regions', data['user_00000'])

        tsv_lines = first['tsv'].decode().splitlines()
        primary = tsv_lines[0].split('\t')
        secondary = tsv_lines[1].split('\t')
        # Ensure all primary header names are unique for datavzrd
        self.assertEqual(len(primary), len(set(primary)))
        self.assertEqual(primary[:3], ['sample_id', 'vcf_sample_id', 'QC_status'])
        self.assertIn('regions-sn-number_of_snps', primary)
        self.assertIn('bcftools_stats_regions', secondary)

        flags_df = pl.read_csv(self.root / 'var_report.qc_flags.tsv', separator='\t')
        self.assertEqual(flags_df['QC_status'].to_list(), ['PASS', 'PASS', 'FAIL'])

        cfg = yaml.safe_load(first['qc_var_datavzrd.yml'])
        self.assertIn('var_qc_metrics', cfg['datasets'])
        self.assertIn('var_qc_flags', cfg['datasets'])
        self.assertIn('var_qc_criteria', cfg['datasets'])
        self.assertIn('var_qc_metrics', cfg['views'])

        # Run without optional regions manifest
        no_reg_root = self.root / 'no_regions'
        no_reg_root.mkdir()
        args_no_reg = make_cohort(no_reg_root, samples=2, include_regions=False)
        self.run_main(args_no_reg)
        primary_no_reg = (no_reg_root / 'var_report.tsv').read_text().splitlines()[0].split('\t')
        self.assertNotIn('regions-sn-number_of_snps', primary_no_reg)


if __name__ == '__main__':
    unittest.main()
