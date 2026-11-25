# Development & Testing

## Running Module Tests (`nf-test`)

All local Nextflow modules under `modules/local/` include `nf-test` suites (`tests/main.nf.test`) covering both standard execution and `-stub` mode. Shared test configuration (`params.modules_testdata_base_path` and single-CPU resource limits) is defined in `nf-test.config` and `conf/nf-test.config`.

Run all local module tests with Singularity:

```bash
nf-test test modules/local/**/tests/main.nf.test --profile singularity
```

Run a single module test:

```bash
nf-test test modules/local/prepareqcreport/tests/main.nf.test --profile singularity
```

## Running Python Unit Tests

The reporting scripts in `bin/` (`summarise_align_qc.py`, `summarise_var_qc.py`, `create_qc_report.py`) have unit test suites under `tests/`:

```bash
pytest tests/
```

## Previewing Documentation Locally

Documentation is built with [MkDocs](https://www.mkdocs.org/) and [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/):

```bash
pip install mkdocs-material
mkdocs serve
```

Then open `http://127.0.0.1:8000` in your browser. To run a strict static build matching CI:

```bash
mkdocs build --strict
```

## Automated GitHub Pages Publishing

The GitHub Actions workflow `.github/workflows/docs.yml` automatically builds and deploys the MkDocs site to GitHub Pages whenever a **GitHub Release** is published (and can also be triggered manually via **Actions → Build and Deploy Docs → Run workflow**).

To enable this in your repository settings:

1. Go to **Settings → Pages** in GitHub.
2. Under **Build and deployment → Source**, select **GitHub Actions**.
