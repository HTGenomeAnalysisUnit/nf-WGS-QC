# Running the Pipeline

## Execution Modes

`nf-WGS-QC` supports running per-sample QC and cohort reporting cleanly in separate invocations:

- **Alignment QC only**: pass `--input_bam`
- **Variant QC only**: pass `--input_vcf`
- **Alignment + Variant QC together**: pass both `--input_bam` and `--input_vcf`
- **Cohort QC Reports**: pass `--make_reports true` (must be run separately; cannot be combined with `--input_bam` or `--input_vcf`)

!!! note "Why report generation is separate"
    Separating `MAKE_REPORTS` allows you to run per-sample QC on individual batches as data arrives, then aggregate one or multiple batches into a single cohort-wide Somalier relatedness/ancestry analysis and unified `datavzrd` report.

---

## Using a Parameters File

For reproducibility, store parameters in a YAML file (`params.yml`) and pass it with `-params-file`:

=== "Alignment & Variant QC"

    ```yaml title="params_qc.yml"
    project_id: "my_wgs_project"
    outdir: "/path/to/QC_results"
    genome_fasta: "/ref/GRCh38.fasta"
    regions_bed: "/ref/callable_regions.bed.gz"

    # Alignment QC
    input_bam: "samples_bam.tsv"
    somalier_vars: "/ref/somalier/sites.hg38.vcf.gz"
    verifybamid_svd_file_prefix: "/ref/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat"

    # Variant QC
    input_vcf: "samples_vcf.tsv"
    af_annotation_zip: "/ref/echtvar/gnomad.v4.echtvar.zip"
    ```

    ```bash
    nextflow run HTGenomeAnalysisUnit/nf-WGS-QC \
        -profile singularity \
        -params-file params_qc.yml
    ```

=== "Cohort QC Reports"

    ```yaml title="params_reports.yml"
    project_id: "my_wgs_project"
    outdir: "/path/to/QC_results"
    make_reports: true

    # Either a single project folder...
    project_qc_folder: "/path/to/QC_results"
    # ...or multi-batch folder lists:
    # alignment_qc_inputs: "align_folders.txt"
    # variant_qc_inputs: "var_folders.txt"

    alignment_qc_criteria: "assets/align_qc_flags.yml"
    variant_qc_criteria: "assets/var_qc_flags.yml"
    somalier_labels: "/ref/somalier/ancestry-labels-1kg.tsv"
    somalier_1kgdata: "/ref/somalier/1kg-somalier"
    ```

    ```bash
    nextflow run HTGenomeAnalysisUnit/nf-WGS-QC \
        -profile singularity \
        -params-file params_reports.yml
    ```

---

## Container & Execution Profiles

| Profile | Description |
| :--- | :--- |
| `singularity` | Runs processes in Singularity containers with `autoMounts = true` (recommended for HPC). |
| `apptainer` | Runs processes in Apptainer containers with `autoMounts = true`. |
| `docker` | Runs processes in Docker containers with `-u $(id -u):$(id -g)`. |
| `conda` / `mamba` | Uses Conda/Mamba environments where supported. **Note**: `QUICK_READS_STATS` and `COMPUTE_CHARR` require containers and will error under `conda`/`mamba`. |
| `test` | Caps process resource limits (`cpus: 1`, `memory: 8.GB`, `time: 1.h`) for lightweight testing. |

### Institutional Profiles

The pipeline automatically loads [nf-core/configs](https://github.com/nf-core/configs) unless `NXF_OFFLINE` is set. You can combine your institutional profile with a container profile:

```bash
nextflow run HTGenomeAnalysisUnit/nf-WGS-QC \
    -profile <institution>,singularity \
    -params-file params_qc.yml
```

### Printing Built-in Help

```bash
nextflow run HTGenomeAnalysisUnit/nf-WGS-QC --help
```
