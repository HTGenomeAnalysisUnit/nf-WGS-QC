"""Regression tests for create_qc_report.py; run with unittest."""

import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest

import polars as pl
import yaml

SCRIPT = Path(__file__).resolve().parents[1] / 'bin' / 'create_qc_report.py'
SPEC = importlib.util.spec_from_file_location('create_qc_report', SCRIPT)
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)


def _sample_align_config(metrics_path, criteria_path):
    return {
        'name': 'QC metrics report',
        'datasets': {
            'qc_metrics': {'path': str(metrics_path), 'separator': '\t', 'headers': 2},
            'qc_flags': {'path': 'align.qc_flags.tsv', 'separator': '\t'},
            'qc_criteria': {'path': str(criteria_path), 'separator': '\t'},
        },
        'views': {
            'qc_metrics': {'dataset': 'qc_metrics', 'desc': 'Alignment QC metrics', 'render-table': {'columns': {}}},
            'qc_flags': {'dataset': 'qc_flags', 'desc': 'Alignment QC flags', 'render-table': {'columns': {}}},
            'qc_criteria': {'dataset': 'qc_criteria', 'desc': 'Alignment QC criteria', 'render-table': {'columns': {}}},
        },
    }


def _sample_var_config(metrics_path, criteria_path):
    return {
        'name': 'QC metrics report',
        'datasets': {
            'var_qc_metrics': {'path': str(metrics_path), 'separator': '\t', 'headers': 2},
            'var_qc_flags': {'path': 'var.qc_flags.tsv', 'separator': '\t'},
            'var_qc_criteria': {'path': str(criteria_path), 'separator': '\t'},
        },
        'views': {
            'var_qc_metrics': {'dataset': 'var_qc_metrics', 'desc': 'Variant QC metrics', 'render-table': {'columns': {}}},
            'var_qc_flags': {'dataset': 'var_qc_flags', 'desc': 'Variant QC flags', 'render-table': {'columns': {}}},
            'var_qc_criteria': {'dataset': 'var_qc_criteria', 'desc': 'Variant QC criteria', 'render-table': {'columns': {}}},
        },
    }


class CreateQCReportTests(unittest.TestCase):
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
            report.main(args)

    def test_combined_align_and_var_report_summary_flags_barplot_and_links(self):
        align_flags = self.write(
            'align.qc_flags.tsv',
            'sample_id\tbam_sample_id\tQC_status\tN_WARN_conditions\tN_FAIL_conditions\tWARN-dup\tWARN-never\tFAIL-low_cov\tFAIL-never_fail\n'
            '001\tB001\tPASS\t0\t0\t0\t0\t0\t0\n'
            '002\tB002\tWARN\t1\t0\t1\t0\t0\t0\n'
            '003\tB003\tFAIL\t0\t1\t0\t0\t1\t0\n'
            '004\tB004\tPASS\t0\t0\t0\t0\t0\t0\n',
        )
        var_flags = self.write(
            'var.qc_flags.tsv',
            'sample_id\tvcf_sample_id\tQC_status\tN_WARN_conditions\tN_FAIL_conditions\tWARN-low_tstv\tFAIL-high_charr\tFAIL-never_var\n'
            '001\tV001\tPASS\t0\t0\t0\t0\t0\n'
            '002\tV002\tPASS\t0\t0\t0\t0\t0\n'
            '003\tV003\tFAIL\t0\t1\t0\t1\t0\n'
            '004\tV004\tFAIL\t0\t1\t0\t1\t0\n',
        )
        align_metrics = self.write('align.tsv', 'sample_id\n')
        var_metrics = self.write('var.tsv', 'sample_id\n')
        align_criteria = self.write('align.qc_criteria.tsv', 'severity\tflag\tdescription\n')
        var_criteria = self.write('var.qc_criteria.tsv', 'severity\tflag\tdescription\n')
        align_yml = self.write('align.yml', yaml.safe_dump(_sample_align_config('old_align.tsv', align_criteria)))
        var_yml = self.write('var.yml', yaml.safe_dump(_sample_var_config('old_var.tsv', var_criteria)))

        out_prefix = self.root / 'proj'
        self.run_main([
            '--project_id', 'MYPROJ',
            '--align_datavzrd_yml', str(align_yml),
            '--align_metrics_tsv', str(align_metrics),
            '--align_qc_flags_tsv', str(align_flags),
            '--var_datavzrd_yml', str(var_yml),
            '--var_metrics_tsv', str(var_metrics),
            '--var_qc_flags_tsv', str(var_flags),
            '--output_prefix', str(out_prefix),
        ])

        summary = pl.read_csv(f'{out_prefix}.project_qc_summary.tsv', separator='\t')
        self.assertEqual(summary['section'].to_list(), ['alignment_qc', 'variants_qc', 'total'])
        self.assertEqual(summary['N_samples'].to_list(), [4, 4, 4])
        self.assertEqual(summary['N_PASS'].to_list(), [2, 2, 1])
        self.assertEqual(summary['N_WARN'].to_list(), [1, 0, 1])
        self.assertEqual(summary['N_FAIL'].to_list(), [1, 2, 2])
        self.assertEqual(summary['failure_rate'].to_list(), [0.25, 0.5, 0.5])

        merged = pl.read_csv(f'{out_prefix}.qc_flags.tsv', separator='\t', schema_overrides={
            'sample_id': pl.String, 'bam_sample_id': pl.String, 'vcf_sample_id': pl.String,
        })
        self.assertEqual(
            merged.columns,
            [
                'sample_id', 'bam_sample_id', 'vcf_sample_id',
                'N_warning_flags', 'N_fail_flags', 'QC_status',
                'ALIGN-WARN-dup', 'ALIGN-FAIL-low_cov', 'VAR-FAIL-high_charr',
            ],
        )
        # Active flag columns are encoded as categorical strings so datavzrd treats them as categories
        self.assertEqual(merged['ALIGN-WARN-dup'].to_list(), ['PASS', 'WARN', 'PASS', 'PASS'])
        self.assertEqual(merged['ALIGN-FAIL-low_cov'].to_list(), ['PASS', 'PASS', 'FAIL', 'PASS'])
        self.assertEqual(merged['VAR-FAIL-high_charr'].to_list(), ['PASS', 'PASS', 'FAIL', 'FAIL'])

        # All 7 configured flags appear in qc_flags_summary.tsv (2 rows each: PASS + FAIL/WARN)
        flags_summary = pl.read_csv(f'{out_prefix}.qc_flags_summary.tsv', separator='\t')
        self.assertEqual(flags_summary.height, 14)
        high_charr_rows = flags_summary.filter(pl.col('flag') == 'VAR-FAIL-high_charr')
        self.assertEqual(high_charr_rows['status'].to_list(), ['PASS', 'FAIL'])
        self.assertEqual(high_charr_rows['count'].to_list(), [2, 2])
        dup_rows = flags_summary.filter(pl.col('flag') == 'ALIGN-WARN-dup')
        self.assertEqual(dup_rows['status'].to_list(), ['PASS', 'WARN'])
        self.assertEqual(dup_rows['count'].to_list(), [3, 1])

        cfg = yaml.safe_load(Path(f'{out_prefix}.datavzrd.yml').read_text())
        self.assertEqual(cfg['name'], 'MYPROJ QC report')
        self.assertEqual(cfg['default-view'], '01_project_qc_summary')
        self.assertEqual(cfg['max-in-memory-rows'], 50000)
        summary_js = Path(f'{out_prefix}.project_qc_summary.js')
        self.assertTrue(summary_js.exists())
        self.assertEqual(cfg['views']['01_project_qc_summary']['render-html']['script-path'], str(summary_js))
        self.assertEqual(
            list(cfg['views'].keys()),
            [
                '01_project_qc_summary',
                '02_qc_flags',
                '03_qc_flags_summary',
                '04_align_qc_metrics',
                '05_var_qc_metrics',
                '06_align_qc_criteria',
                '07_var_qc_criteria',
            ],
        )
        spec = json.loads(cfg['views']['03_qc_flags_summary']['render-plot']['spec'])
        self.assertEqual([p['title'] for p in spec['vconcat']], ['FAIL flags', 'WARN flags'])
        self.assertEqual(spec['vconcat'][0]['encoding']['color']['scale']['domain'], ['PASS', 'FAIL'])
        self.assertEqual(spec['vconcat'][1]['encoding']['color']['scale']['domain'], ['PASS', 'WARN'])
        self.assertEqual(
            cfg['views']['02_qc_flags']['render-table']['columns']['ALIGN-WARN-dup']['plot']['heatmap']['domain'],
            ['PASS', 'WARN'],
        )
        self.assertEqual(
            cfg['views']['02_qc_flags']['render-table']['columns']['VAR-FAIL-high_charr']['plot']['heatmap']['domain'],
            ['PASS', 'FAIL'],
        )
        self.assertEqual(
            cfg['datasets']['02_qc_flags']['links']['04_align_qc_metrics'],
            {'column': 'sample_id', 'table-row': '04_align_qc_metrics/sample_id', 'optional': True},
        )

    def test_custom_max_rows_per_page_and_align_only_mode(self):
        align_flags = self.write(
            'align.qc_flags.tsv',
            'sample_id\tbam_sample_id\tQC_status\tN_WARN_conditions\tN_FAIL_conditions\tFAIL-low_cov\n'
            'S1\tB1\tFAIL\t0\t1\t1\n'
            'S2\tB2\tPASS\t0\t0\t0\n',
        )
        align_yml = self.write('align.yml', yaml.safe_dump(_sample_align_config('a.tsv', 'ac.tsv')))

        align_only_prefix = self.root / 'align_only'
        self.run_main([
            '--max_rows_per_page', '12345',
            '--align_datavzrd_yml', str(align_yml),
            '--align_qc_flags_tsv', str(align_flags),
            '--output_prefix', str(align_only_prefix),
        ])
        cfg_a = yaml.safe_load(Path(f'{align_only_prefix}.datavzrd.yml').read_text())
        self.assertEqual(cfg_a['max-in-memory-rows'], 12345)
        self.assertEqual(
            list(cfg_a['views'].keys()),
            ['01_project_qc_summary', '02_qc_flags', '03_qc_flags_summary', '04_align_qc_metrics', '06_align_qc_criteria'],
        )
        spec_a = json.loads(cfg_a['views']['03_qc_flags_summary']['render-plot']['spec'])
        # Only FAIL flags configured -> single FAIL panel
        self.assertEqual([p['title'] for p in spec_a['vconcat']], ['FAIL flags'])


if __name__ == '__main__':
    unittest.main()
