    # Save a summary table reporting the overall number of PASS, WARN and FAIL samples
    summary_df = df['QC_status'].value_counts().sort('QC_status')
    summary_df.write_csv(f"{args.output_prefix}.qc_summary.tsv", separator='\t')
    log(f"QC summary file saved to {args.output_prefix}.qc_summary.tsv")