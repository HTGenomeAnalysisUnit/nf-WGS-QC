# QC Variants

Generate QC report from genomic variants VCF

## Select files by range date

```bash
find -newerct "1 Aug 2013" ! -newerct "1 Sep 2013" -ls
```

Female = 0
Male = 1

## Implementation

1. QC metrics are automatically generated as new BAM / VCF files appears in watched folder
2. Make tranches across date intervals using `make_tranche_time.py`
  
    ```bash
    make_tranches_time.py --manifest_file bam_manifest.tsv --start_date 20221205
    ```

3. Use the `MAKE_REPORTS` entrypoint to generate reports given the tranches data
4. Finally, merge VCF files from a tranche and generate variants report per tranche by running the Rmd report template
