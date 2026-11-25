# Variant QC Workflow (`VARIANTS_QC`)

The `VARIANTS_QC` workflow (`workflows/qc_variants.nf`) runs when `--input_vcf` is provided. It handles standard VCFs, gVCFs, and BCFs and publishes per-sample metrics plus TSV manifests for `MAKE_REPORTS`.

## Pipeline Steps

1. **Input validation & indexing (`BCFTOOLS_INDEX`)**
   - Reads and validates `--input_vcf` against `assets/schema_input_vcf.json`.
   - Detects file type (`gvcf` vs `vcf`) and existing `.tbi`/`.csi` indexes; creates missing indexes via `bcftools index`.
2. **Variant filtering (`BCFTOOLS_VIEW`)**
   - Filters input variants according to `--variant_stats_filter` (default `.,PASS`) and normalizes representation for downstream stats and contamination estimation.
3. **VCF sample manifest (`VCF_MANIFEST`)**
   - Extracts internal VCF sample IDs, file dates, and paths into `input_vcf_manifest.tsv` and `sample_to_user_map.tsv`, mapping each user sample ID to its internal `vcf_sample_id`.
4. **Contamination estimation (`VCF_CONTAMINATION` subworkflow)**
   - When `--annotate_af true` (default):
     1. Normalizes variants with `BCFTOOLS_NORM` against `--genome_fasta`.
     2. Annotates population allele frequencies with `ECHTVAR_ANNO` using `--af_annotation_zip`.
     3. Indexes the annotated VCF (`BCFTOOLS_INDEX`).
   - Runs `COMPUTE_CHARR` (`sceVCF`) using the configured `--charr_af_field` (default `gnomad_af`) to estimate variant-level contamination (`CHARR`) and inconsistent heterozygous allele-balance rates.
   - Requires a container runtime (Singularity, Apptainer, or Docker).
5. **Variant statistics (`BCFTOOLS_STATS`, `BCFTOOLS_IDXSTATS`)**
   - `BCFTOOLS_STATS` runs on all filtered variants to collect SNV/indel counts, Ts/Tv ratios, QUAL distributions, and genotype depth (`DP`) distributions.
   - When a regions BED is provided (via the samplesheet or `--regions_bed`), `BCFTOOLS_STATS_WITH_REGIONS` additionally computes stats restricted to those regions.
   - `BCFTOOLS_IDXSTATS` computes per-chromosome variant counts from the VCF index.
