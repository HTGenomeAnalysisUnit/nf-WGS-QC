# Alignment QC Workflow (`ALIGNMENTS_QC`)

The `ALIGNMENTS_QC` workflow (`workflows/qc_alignments.nf`) runs when `--input_bam` is provided. It processes BAM and CRAM files in parallel and publishes both per-sample outputs and Nextflow-indexed TSV manifests consumed by `MAKE_REPORTS`.

## Pipeline Steps

1. **Input validation & indexing (`SAMTOOLS_INDEX`)**
   - Reads and validates `--input_bam` against `assets/schema_input_bam.json`.
   - Automatically detects existing `.bai` / `.crai` indexes; unindexed files are indexed with `samtools index`.
2. **BAM sample manifest (`BAM_MANIFEST`)**
   - Extracts internal BAM/CRAM sample IDs (`@RG SM`), file modification dates, and paths into `input_bam_manifest.tsv` and `sample_to_user_map.tsv`.
   - Joins the extracted `bam_sample_id` onto every downstream output's metadata so cohort reporting can reconcile user IDs with internal read-group sample IDs.
3. **Coverage profiling (`MOSDEPTH`)**
   - Computes global coverage distribution (`global.dist.txt`) and per-contig summary (`mosdepth.summary.txt`).
   - When a regions BED is provided (via the samplesheet or `--regions_bed`), computes regional coverage and base fractions above `--mosdepth_thresholds`.
   - Produces quantized coverage beds using `--mosdepth_quantize_bins` and `--mosdepth_quantize_labels`.
   - Optionally computes per-base coverage when `--mosdepth_per_base_coverage true`.
4. **Alignment statistics (`SAMTOOLS_STATS`, `SAMTOOLS_IDXSTATS`)**
   - `samtools stats` computes read counts, mapping rate, error rate, insert size, and base quality distributions using `--genome_fasta`.
   - `samtools idxstats` computes per-chromosome mapped/unmapped read counts for autosomal uniformity and coverage checks.
5. **Per-flowcell / lane read statistics (`QUICK_READS_STATS`)**
   - Runs `qrs` on up to `--qrs_subset_reads` reads (default `5,000,000`) to summarize sequencing run, flowcell, lane, and read-level quality metrics (`per_sample_run.tsv`, `rundetails.tsv`, `read_stats.tsv`).
   - Requires a container runtime (Singularity, Apptainer, or Docker).
6. **Contamination estimation (`VERIFYBAMID_VERIFYBAMID2`)**
   - Estimates DNA sample contamination (`FREEMIX`) and genetic PCs against the provided SVD panel (`--verifybamid_svd_file_prefix`).
7. **Genetic fingerprinting (`SOMALIER_EXTRACT`)**
   - Extracts informative genotypes at `--somalier_vars` sites into compact `.somalier` binary profiles for cohort-level sex, relatedness, and ancestry checks in `MAKE_REPORTS`.
