# QC Variants

Generate QC report from genomic variants VCF

## Select files by range date

```bash
find -newerct "1 Aug 2013" ! -newerct "1 Sep 2013" -ls
```

Female = 0
Male = 1

New implementation

1. Generate all metrics files as new BAM/VCF files are loaded using a watchPath channel. This output files in a general QC folder under qc_variants and qc_alignments.
2. The above emits a 2 col list of samples ids in the general QC folder that we update in the app/data/project dir. The app then read this to give real time number of downloaded samples.
3. Make a new make_report workflow that read files from the general QC folder, partition them in time intervals and then runs somalier ancetry/profiles, multiQC and finally the R script to generate reports, summary qc table and RDS objects for the app. This is run manually once a week.
4. The above workflow updates the app loading the new RDS

- Can we update the app using rsconnect but loading only the RDS files in app_folder/data?
- App now reports N available samples (real-time update) and N processed samples (update per week). Stats and plots refers only to processed samples using the RDS objects.
- App must read FAIL/PASS status in summary QC and provide Infobbxes based on this
- App must read summary QC and provide per tranche stats on N samples and % failed samples.
