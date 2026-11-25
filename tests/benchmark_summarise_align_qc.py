"""Generate a cohort and time report generation (input generation is excluded).

python tests/benchmark_summarise_align_qc.py --samples 1000 --workers 1
Use --script /path/to/another/version.py to compare implementations.
Peak RSS covers the entire process, including fixture creation and validation.
"""

import argparse
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import resource
import tempfile
import time

from test_summarise_align_qc import make_cohort, SCRIPT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=1000)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--cycles', type=int, default=150)
    parser.add_argument('--max-depth', type=int, default=200)
    parser.add_argument('--script', type=Path, default=SCRIPT)
    options = parser.parse_args()
    if min(options.samples, options.workers, options.cycles, options.max_depth) < 1:
        parser.error('All counts must be positive')
    spec = importlib.util.spec_from_file_location('benchmark_target', options.script)
    target = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(target)
    target.log = lambda *args, **kwargs: None
    with tempfile.TemporaryDirectory(prefix='align-qc-benchmark-') as directory:
        args = make_cohort(Path(directory), options.samples, options.cycles, options.max_depth)
        if options.workers != 1:
            args += ['--workers', str(options.workers)]
        started = time.perf_counter()
        # Do not retain one Python string per logged metric file.
        with contextlib.redirect_stdout(io.StringIO()):
            target.main(args)
        elapsed = time.perf_counter() - started
        result = json.loads((Path(directory) / 'report.json').read_text())
        assert len(result) == options.samples
    print(json.dumps({
        'samples': options.samples, 'metric_files': options.samples * 8,
        'somalier_pairs': options.samples * (options.samples - 1) // 2,
        'workers': options.workers, 'seconds': round(elapsed, 3),
        'peak_rss_mib': round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
    }, indent=2))


if __name__ == '__main__':
    main()
