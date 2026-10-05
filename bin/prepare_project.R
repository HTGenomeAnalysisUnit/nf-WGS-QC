#!/usr/bin/env Rscript

library(rsconnect)
library(optparse)
library(rlang)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# Constants --------------

min_cov <- 10
min_qual <- 20
canonical_chrs <- paste("chr",c(1:22,"X","Y","M"), sep="")
timestamp <- format(Sys.time(), "%Y/%m/%d %H:%M:%S")
file_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Helper functions -----------------
`%nin%` <- Negate('%in%')

save_tsv <- function(df, outfile) {
  write.table(df, file = outfile, sep="\t", quote=F, row.names = F)
}

evaluate_filters <- function(df, filters_file, tag) {
  if (!is.na(filters_file)) {
    df$tmp_col_reason <- "PASS"
    #df$QC_pass_status <- "PASS"
    for (filter_string in scan(filters_file, what="", sep="\n", comment.char="#")) {
      df <- df %>% 
        mutate(
          QC_pass_status = ifelse(
            eval_tidy(parse_expr(filter_string)), 
            tag, 
            QC_pass_status
          ),
          tmp_col_reason = ifelse(
            eval_tidy(parse_expr(filter_string)), 
            paste(tmp_col_reason,filter_string,sep=";"), 
            tmp_col_reason
          )
        ) 
    }
    df <- df %>% mutate(
        tmp_col_reason=ifelse(
          tmp_col_reason=="PASS",
          "PASS",
          gsub("PASS;","",tmp_col_reason)
        )
      )
    reason_col_name <- paste0(tag, "_reason")
    df <- df %>% rename( !!reason_col_name := tmp_col_reason )
  }
  return(df)
}

merge_tables <- function(tables, by="sample") {
  res <- tables[1]
  for (t in tables[-1]) {
    res <- merge(res, t, by=by, all.x=T)
  }
  return(res)
}

make_min_max_value <- function(x, scaling, precision) {
  s <- paste0(
    round(min(x, na.rm = T)/scaling, digits = precision), "-", 
    round(max(x, na.rm = T)/scaling, digits = precision)
  )
  return(s)
}

checkNumbers <- function(values) {
  if (length(grep("^[0-9.]+$", values, perl=T)) == length(values)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Variants QC functions ---------
readVCFstat <- function(stat_file, sep="\t", fileID=FALSE) {
  outdata <- list(
    SN = NULL,
    TSTV = NULL,
    SiS = NULL,
    AF = NULL,
    QUAL = NULL,
    IDD = NULL,
    ST = NULL,
    DP = NULL,
    PSC = NULL,
    PSI = NULL,
    HWE = NULL 
    )
  
  headers <- list (
    SN = c("SN","id","key","value"),
    TSTV = c("TSTV", "id", "ts", "tv", "ts_tv", "ts_1st_ALT","tv_1st_ALT", "ts_tv_1st_ALT"),
    SiS = c("SiS", "id", "AC", "n_SNPs", "n_transitions", "n_transversions", "n_indels", "repeat_consistent", "repeat_inconsistent", "not_applicable"),
    AF = c("AF", "id", "Allele_freq", "n_SNPs", "n_transitions", "n_transversions", "n_indels", "repeat_consistent", "repeat_inconsistent", "not_applicable"),
    QUAL = c("QUAL", "id", "quality", "n_SNPs", "n_transitions_1stALT", "n_transversions_1stALT", "n_indels"),
    IDD = c("IDD", "id", "length", "n_sites", "n_genotypes", "mean_VAF"),
    ST = c("ST", "id", "type", "count"),
    DP = c("DP", "id", "bin", "n_genotypes", "pct_genotypes", "n_sites", "pct_sites"),
    PSC = c("PSC", "id", "sample", "nRefHom", "nNonRefHom", "nHets", "nTransitions", "nTransversions", "nIndels", "avg_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissing"),
    PSI = c("PSI", "id", "sample", "in_frame", "out_frame", "not_applicable", "out_in.out_ratio", "nInsHets", "nDelHets", "nInsAltHoms", "nDelAltHoms"),
    HWE = c("HWE", "id", "1st_ALT_AF", "n_observations", "25th_percentile", "median", "75th_percentile")
  )
  
  
  con <- file(stat_file, "r")
  lines <- c()
  reading_tag <- "START"
  while(TRUE) {
    line = readLines(con, 1)
    if(length(line) == 0) {
      break
    } else {
      if (startsWith(line, "#")) {next}
      line <- unlist(strsplit(line, sep))
      tag <- line[1]
      if (reading_tag != tag) {
        reading_tag <- tag  
      }
      switch(tag,
             SN = {outdata$SN <- rbind(outdata$SN, line)},
             TSTV = {outdata$TSTV <- rbind(outdata$TSTV, line)},
             SiS = {outdata$SiS <- rbind(outdata$SiS, line)},
             AF = {outdata$AF <- rbind(outdata$AF, line)},
             QUAL = {outdata$QUAL <- rbind(outdata$QUAL, line)},
             IDD = {outdata$IDD <- rbind(outdata$IDD, line)},
             ST = {outdata$ST <- rbind(outdata$ST, line)},
             DP = {outdata$DP <- rbind(outdata$DP, line)},
             PSC = {outdata$PSC <- rbind(outdata$PSC, line)},
             PSI = {outdata$PSI <- rbind(outdata$PSI, line)},
             HWE = {outdata$HWE <- rbind(outdata$HWE, line)} )
    }
  }
  
  for (n in names(outdata)) {
    if (!is.null(outdata[[n]])) {
      outdata[[n]] <- as.data.frame(outdata[[n]], stringsAsFactors=F)
      colnames(outdata[[n]]) <- headers[[n]]
      
      
      outdata[[n]] <- as.data.frame(outdata[[n]] %>% mutate_if(checkNumbers,as.numeric))
    }
  }
  
  if (fileID != FALSE) {
    for (n in names(outdata)) {
      if (!is.null(outdata[[n]])) {outdata[[n]]$sample <- fileID}
    }
  }
  close(con)
  return(outdata) 
}

read_contig_vars <- function(path, contigs) {
  message("### Processing contig count data ###")
  message("### Path: ", path)
  res <- list(by_sample=list(), all_samples=NULL)
  
  files <- list.files(path, pattern = ".PASS.varsPerContig.tsv")
  for (f in files) {
    message("Processing ", f)
    sample_name <- gsub("\\.PASS\\.varsPerContig\\.tsv","", f)
    dt <- fread(paste0(path, "/", f), header=F, col.names=c("contig","size","n_PASS_vars"))
    if (nrow(dt) == 0) {
      message("No data for ", sample_name, "; setting to zero")
      dt <- data.table(contig=contigs, size=0, n_PASS_vars=0, vars_per_Mb=0, sample=sample_name)
    } else {
      idx <- which(dt$contig %in% contigs)
      dt <- dt[idx]
      dt$vars_per_Mb <- dt$n_PASS_vars / (dt$size / 1000000)
      # In case some of the requested contigs is missing in stats file we add them with zero counts
      missing_contigs <- contigs[!contigs %in% dt$contig]
      if (length(missing_contigs) > 0) {
        message("Missing contigs for ", sample_name, ": ", paste(missing_contigs, collapse=","))
        missing_dt <- data.table(contig=missing_contigs, size=0, n_PASS_vars=0, vars_per_Mb=0)
        dt <- rbind(dt, missing_dt)
      }
      # Set the sample name
      dt$sample <- sample_name
    }
    res$by_sample[[sample_name]] <- dt
    res$all_samples <- rbind(res$all_samples, dt)
  }
  return(res)
}

read_stats_vars <- function(path, minDP, minQUAL) {
  message("### Processing vars stats data ###")
  message("### Path: ", path)
  
  res <- list(by_sample=list(), DP=NULL, QUAL=NULL, summary=NULL)
  
  files <- list.files(path, pattern = "vars.stats")
  for (f in files) {
    message("Processing ", f)
    sample_name <- gsub("\\.vars\\.stats","", f)
    sdata <- readVCFstat(paste0(path, "/", f), fileID = sample_name)
    res$by_sample[[sample_name]] <- sdata

    # Do not process sample data if any of the stats sections are empty
    if (any(sapply(sdata, is.null))) {
      message("Skipping sample ", sample_name, " due to missing stats data")
      next
    }
    #Compile a summary tab across all samples
    genos_above_mindp <- sdata$DP %>% 
      arrange(desc(as.numeric(bin))) %>% 
      mutate(cumsum = cumsum(n_genotypes), pct_cumsum = cumsum/sum(n_genotypes)) %>% 
      filter(bin == minDP) %>% select(sample, pct_cumsum) %>% 
      rename_(.dots=setNames("pct_cumsum", paste0("pct_genos_dp_above_",minDP)))
    genos_above_minqual <- sdata$QUAL %>% 
      arrange(desc(as.numeric(quality))) %>% 
      mutate(
        cumsum_snp = cumsum(n_SNPs), 
        pct_snp = cumsum_snp/sum(n_SNPs),
        cumsum_indels = cumsum(n_indels), 
        pct_indels = cumsum_indels/sum(n_indels)
      ) %>% 
      filter(quality == minQUAL) %>% select(sample, pct_snp, pct_indels) %>% 
      rename_(.dots=setNames(c("pct_snp","pct_indels"), c(paste0("pct_snp_qual_above_",minQUAL),paste0("pct_indels_qual_above_",minQUAL))))
    vars_counts <- sdata$PSC %>% select(-PSC, -id) %>%
      mutate(nSNPs = nNonRefHom + nHets, Ts_Tv = nTransitions/nTransversions) %>%
      select(sample, nSNPs, nIndels, nNonRefHom, nHets, Ts_Tv, nMissing, avg_depth) %>%
      rename(nSNPsAltHom = nNonRefHom, nSNPsHets = nHets)
    indels_counts <- sdata$PSI %>% 
      select(sample, nInsHets, nDelHets, nInsAltHoms, nDelAltHoms)
    summary_dt <- merge_tables(
      list(vars_counts,
      indels_counts,
      genos_above_mindp,
      genos_above_minqual),
      by = "sample"
    )
    
    res$summary <- rbind(res$summary,summary_dt)
    
    if (is.null(res$DP)) {
      res$DP <- sdata$DP %>% select(bin, n_genotypes) 
    } else {
    res$DP <- merge(
      res$DP,
      sdata$DP %>% select(bin, n_genotypes) %>% rename(new_genos = n_genotypes),
      by = "bin", all=T) %>%
      replace_na(list(n_genotypes = 0, new_genos = 0)) %>%
      mutate(n_genotypes = n_genotypes + new_genos) %>%
      select(-new_genos)
    }
    
    if (is.null(res$QUAL)) {
      res$QUAL <- sdata$QUAL %>% select(quality, n_SNPs, n_indels) 
    } else {
      res$QUAL <- merge(
        res$QUAL,
        sdata$QUAL %>% select(quality, n_SNPs, n_indels) %>% 
          rename(new_snps = n_SNPs, new_indels = n_indels),
        by = "quality", all=T) %>%
        replace_na(list(n_SNPs = 0, new_snps = 0, n_indels = 0, new_indels = 0)) %>%
        mutate(
          n_SNPs = n_SNPs + new_snps,
          n_indels = n_indels + new_indels) %>%
        select(-new_snps,-new_indels)
    }
  }
  return(res)
}

read_contamination_data <- function(path) {  
  message("### Processing contamination data ###")
  message("### Path: ", path)
  
  cont_res <- NULL
  files <- list.files(path, pattern = "cont.tsv")
  for (f in files) {
    message("Processing ", f)
    dt <- fread(paste0(path, "/", f), header=T, check.names=T)
    setnames(dt, "X.SAMPLE", "sample")
    cont_res <- rbind(cont_res, dt)
  }
  return(cont_res)
}

read_vars_QC <- function(project_id, data_path, minDP, minQUAL, by_tranche=T, previous_rds=NA, force_update=F) {
  stats <- list()
  vars_per_contig <- list()
  summary_qc <- NULL
  contamination <- list()
  old_rds <- NULL
  
  # If we need to update a previous RDS read it and clean previous applied filter columns
  if (!is.na(previous_rds)) {
    old_rds <- readRDS(previous_rds)
  }

  tranches <- list()
  if (by_tranche) {
    tranche_name <- list.dirs(data_path, recursive = F,full.names = F)
    tranche_name <- tranche_name[ grepl("tranche", tranche_name) ]
    for (t in tranche_name) {
      tranche_path <- paste0(data_path,"/",t,"/qc_variants")
      tranches[[t]] <- tranche_path
    }
  } else {
    tranches[[project_id]] <- paste0(data_path,"/qc_variants")
  }

  # Load data from a previous RDS is present
  if (!is.null(old_rds)) {
    for (t in names(tranches)) {
      # if we want to force update then we need to remove new tranches that are aldready present
      if (force_update) {
        if (t %in% names(old_rds$coverage)) {
          old_rds$stats[[t]] <- NULL
          old_rds$vars_per_contig[[t]] <- NULL
          old_rds$contamination[[t]] <- NULL
          old_rds$summary_qc %>% filter(group != t)
        }
      } else {
        # If not forcing update than tranches already in the previous RDS will be skipped
        if (t %in% names(old_rds$coverage)) {
          tranches[[t]] <- NULL
        }
      }
    }
    # Store the relevant summary data so we can update the tables
    summary_qc <- old_rds$summary_qc %>% select(-QC_pass_status,-fail_reason)
    stats <- old_rds$stats
    vars_per_contig <- old_rds$vars_per_contig
    contamination <- old_rds$contamination
  } 

  for (t in names(tranches)) {
    tranche_path <- tranches[[t]]
    var_stats_path <- paste0(tranche_path,"/variant_stats")
    contamination_data_path <- paste0(tranche_path,"/contamination_stats")
    if (!dir.exists(var_stats_path) | !dir.exists(contamination_data_path)) {
      message("No variant QC data data found for ", t)
      next
    }
    stats[[t]] <- read_stats_vars(var_stats_path, minDP, minQUAL)
    vars_per_contig[[t]] <- read_contig_vars(var_stats_path, canonical_chrs)
    contamination[[t]] <- read_contamination_data(contamination_data_path)
    no_vars_contigs <- as.data.table(vars_per_contig[[t]]$all_samples %>% 
      filter(n_PASS_vars == 0) %>% group_by(sample) %>% 
      summarise(contigs_with_zero_vars = paste(contig, collapse=",")) %>%
      ungroup())

    list_of_dt <- list(
      stats[[t]]$summary,
      no_vars_contigs,
      contamination[[t]][,.(sample, CHARR, HETEROZYGOSITY_RATE, INCONSISTENT_AB_HET_RATE)]
    )
    merged_dt <- merge_tables(list_of_dt)
    
    #Compute charr mean / sd / deviation across the tranche
    charr_stats <- list(mean=mean(merged_dt$CHARR), sd=sd(merged_dt$CHARR))
    het_rate_stats <- list(mean=mean(merged_dt$HETEROZYGOSITY_RATE), sd=sd(merged_dt$HETEROZYGOSITY_RATE))
    merged_dt$charr_sd_from_mean <- round((merged_dt$CHARR - charr_stats$mean) / charr_stats$sd, 3)
    merged_dt$hetrate_sd_from_mean <- round((merged_dt$HETEROZYGOSITY_RATE - het_rate_stats$mean) / het_rate_stats$sd, 3)
    
    #Add the group tag and rbind
    merged_dt$group <- t
    summary_qc <- rbind(summary_qc, merged_dt)
  }
  message("### VARS QC PROCESSING COMPLETE ###")
  
  return(list(stats=stats, vars_per_contig=vars_per_contig, summary_qc=summary_qc, contamination=contamination))
}

# Align QC functions ----------
flagstat_elements <- c(
"total",
"primary",
"secondary",
"supplementary",
"duplicates",
"primary.duplicates",
"mapped",
"primary.mapped",
"paired",
"read1",
"read2",
"properly.paired",
"itself.and.mate.mapped",
"singletons",
"mate.on.different.chr",
"mate.on.different.chr.mapQ5"
)

read_mosdepth <- function(path, canonical_chrs, group_name, minCov) {
  message("### Processing coverage data ###")
  message("### Path: ", path)
  mincov_dt <- NULL
  cumcov_dt <- NULL
  tot_bases_dt <- NULL

  files <- list.files(path, pattern="totbases")
  for (f in files) {
    message("Processing ", f)
    dt <- fread(paste0(path,"/",f), header=T) 
    dt$mean_cov_bases <- round(dt$N_bases / dt$Genome_size,2)
    tot_bases_dt <- rbind(tot_bases_dt, dt)  
  }

  files <- list.files(path, pattern="global.dist.txt")
  for (f in files) {
    message("Processing ", f)
    sample_name <- gsub("\\.mosdepth\\.global\\.dist\\.txt","",f)
    dt <- fread(cmd = paste0("grep total ", path, "/", f), header=F, 
                col.names=c("min_cov","value"), drop=c(1))
    dt[, sample := sample_name ]
    cumcov_dt <- rbind(cumcov_dt, dt)
    if (max(dt$min_cov) < minCov) {
      message("WARNING: minCov is higher than the max coverage in the file: ", f)
      rbind(dt, list(minCov, 0, sample_name))
    }
    mc <- dt[min_cov == minCov][, !"min_cov"]
    setnames(mc, "value", paste0("pct_bases_min_",minCov))
    mincov_dt <- rbind(mincov_dt, mc)
  }
  cumcov_dt$group <- group_name
  
  meancov_dt <- NULL
  summary_path <- paste0(path, "/summary")
  message("### Path: ", summary_path)
  files <- list.files(summary_path, pattern="*.summary.txt")
  for (f in files) {
    message("Processing ", f)
    sample_name <- gsub("\\.mosdepth\\.summary\\.txt","",f)
    dt <- fread(paste0(summary_path, "/", f), header=T, 
                select=c("chrom","mean","min","max"))
    setnames(dt, old = c("mean","min","max"), new = c("mean_cov","min_cov","max_cov"))
    contigs_idx <- c(which(dt$chrom %in% canonical_chrs), which(dt$chrom == "total"), which(dt$chrom == "autosomes"))
    dt <- dt[contigs_idx]
    dt <- dt[, sample := sample_name ]
    meancov_dt <- rbind(meancov_dt, dt)
  }

  autosomes_cov <- meancov_dt[chrom == "autosomes"]
  autosomes_cov[,chrom := NULL]
  setnames(autosomes_cov, old = c("mean_cov","min_cov","max_cov"), new = c("mean_cov_autosomes","min_cov_autosomes","max_cov_autosomes"))

  meancov_dt$group <- group_name
  total_cov <- meancov_dt[chrom == "total"]
  total_cov[, chrom := NULL]
  
  summary_cov <- merge(total_cov, autosomes_cov, by="sample", all.x=T)
  summary_cov <- merge(summary_cov, mincov_dt, by="sample", all.x=T)
  summary_cov <- merge(summary_cov, tot_bases_dt, by="sample", all.x=T)
  return(list(cumcov=cumcov_dt, meancov=meancov_dt,summarycov=summary_cov))
}

read_samtools_idxstats <- function(path, canonical_chrs, group_name, decoy_string="decoy") {
  res_dt <- NULL
  message("### Processing samtools idxstats data ###")
  message("### Path: ", path)
  files <- list.files(path, pattern="\\.idxstat")
  for (f in files) {
    message("Processing ", f)
    sample_name <- gsub("\\.idxstat","",f)
    dt <- fread(paste0(path, "/", f), header=F, 
                col.names=c("contig","size","n_mapped","n_unmapped"))
    decoys <- grepl(decoy_string,dt$contig)
    decoys_size <- sum(dt$size[decoys])
    decoys_n_mapped <- sum(dt$n_mapped[decoys])
    dt <- dt[!decoys]
    addcontigs <- dt$contig %nin% canonical_chrs
    addcontigs_size <- sum(dt$size[addcontigs])
    addcontigs_n_mapped <- sum(dt$n_mapped[addcontigs])
    dt <- dt[!addcontigs]
    dt <- rbind(dt, 
                list("decoys", decoys_size, decoys_n_mapped, 0),
                list("additional_contigs", addcontigs_size, addcontigs_n_mapped, 0))
    dt[,`:=` (n_mapped_per_mb = n_mapped / (size / 1e6), sample = sample_name) ]
    res_dt <- rbind(res_dt, dt)
  }
  res_dt$group <- group_name
  return(res_dt)
}

read_samtools_flagstats <- function(path, group_name) {
  res_dt <- NULL
  message("### Processing samtools flagstats data ###")
  message("### Path: ", path)
  files <- list.files(path, pattern="\\.flagstat")
  for (f in files) {
    message("Processing ", f)
    sample_name <- gsub("\\.flagstat","",f)
    raw <- read.csv(paste0(path, "/", f), header=F, sep="+")
    values <- raw$V1
    values_perc <- values/values[1]
    dt <- data.frame(metric=c(flagstat_elements,paste0("pct_", flagstat_elements)), N=c(values, values_perc), sample=sample_name)
    res_dt <- rbind(res_dt, dt)
  }
  res_dt <- res_dt %>% pivot_wider(names_from="metric", values_from="N")
  res_dt$group <- group_name
  return(as.data.table(res_dt))
}

read_readstats <- function(path, group_name) {
  res_dt <- NULL
  message("### Processing readstats data ###")
  message("### Path: ", path)
  files <- list.files(path, pattern="\\.per_sample_run\\.tsv")
  for (f in files) {
    message("Processing ", f)
    dt <- fread(paste0(path, "/", f), header=T)
    setnames(dt,"SAMPLE","sample")
    res_dt <- rbind(res_dt, dt)
  }
  res_dt$group <- group_name
  return(res_dt)
}

read_samples_per_run <- function(path, group_name) {
  dt <- NULL
  files <- list.files(path, pattern="\\.rundetails\\.tsv")
  for (f in files) {
    tmp_dt <- fread(paste0(path, "/", f), header=T)
    dt <- rbind(dt, tmp_dt)
  }
  sxrun <- dt[, .(samples_per_run = length(unique(SAMPLE)), id=paste(INSTRUMENT, RUN, sep="_")), by = .(INSTRUMENT, RUN)]
  sxflow <- dt[, .(samples_per_flowcell = length(unique(SAMPLE)), id=paste(INSTRUMENT, RUN, FLOWCELL, sep="_")), by = .(INSTRUMENT, RUN, FLOWCELL)]
  sxlane <- dt[, .(samples_per_lane = length(unique(SAMPLE)), id=paste(INSTRUMENT, RUN, FLOWCELL, LANE, sep="_")), by = .(INSTRUMENT, RUN, FLOWCELL, LANE)]
  sxrun <- sxrun[, 3:4] %>% pivot_longer(samples_per_run, names_to = "metric", values_to = "n")
  sxflow <- sxflow[, 4:5] %>% pivot_longer(samples_per_flowcell, names_to = "metric", values_to = "n")
  sxlane <- sxlane[, 5:6] %>% pivot_longer(samples_per_lane, names_to = "metric", values_to = "n")

  res_dt <- rbind(sxrun,sxflow,sxlane)
  res_dt$group <- group_name
  return(as.data.table(res_dt))
}

read_runs_per_sample <- function(path, group_name) {
  dt <- NULL
  files <- list.files(path, pattern="\\.rundetails\\.tsv")
  for (f in files) {
    tmp_dt <- fread(paste0(path, "/", f), header=T)
    dt <- rbind(dt, tmp_dt)
  }
  runxs <- unique(dt[, .(SAMPLE, INSTRUMENT, RUN)])[, .(runs_per_sample = .N), by = .(SAMPLE)]
  flowxs <- unique(dt[, .(SAMPLE, INSTRUMENT, RUN, FLOWCELL)])[, .(flowcells_per_sample = .N), by = .(SAMPLE)]
  lanexs <- dt[, .(lanes_per_sample = .N), by = .(SAMPLE)]
  runxs <- runxs %>% pivot_longer(runs_per_sample, names_to = "metric", values_to = "n")
  flowxs <- flowxs %>% pivot_longer(flowcells_per_sample, names_to = "metric", values_to = "n")
  lanexs <- lanexs %>% pivot_longer(lanes_per_sample, names_to = "metric", values_to = "n")
  
  res_dt <- rbind(runxs,flowxs,lanexs)
  res_dt$group <- group_name
  res_dt <- as.data.table(res_dt)
  setnames(res_dt, "SAMPLE", "sample")
  return(res_dt)
}

read_somalier <- function(path, group_name) {
  message("### Processing somalier data ###")
  message("### Path: ", path)
  message("Processing ", paste0(path, "/",group_name, ".pairs.tsv"))
  pairs_dt <- fread(paste0(path, "/",group_name, ".pairs.tsv"), header=T)
  setnames(pairs_dt, "#sample_a", "sample_a")
  pairs_dt$group <- group_name
  samples_dt <- fread(paste0(path, "/", group_name, ".samples.tsv"), header=T, 
                      select=c("sample_id","original_pedigree_sex","depth_mean","X_depth_mean","Y_depth_mean","X_het","n_hom_ref","n_het","n_hom_alt"))
  samples_dt$scaled_X_depth <- (samples_dt$X_depth_mean / samples_dt$depth_mean) * 2
  samples_dt$group <- group_name
  setnames(samples_dt, "sample_id", "sample")
  message("Processing ", paste0(path, "/", group_name, ".somalier-ancestry.tsv"))
  ancestry_dt <- fread(paste0(path,"/", group_name, ".somalier-ancestry.tsv"), header=T)
  ancestry_dt$group <- group_name
  ancestry_dt$given_ancestry[ancestry_dt$given_ancestry == ""] <- "External_sample"
  setnames(ancestry_dt, "#sample_id", "sample")
  
  message("Find likely identical pairs")
  related_pairs <- pairs_dt[relatedness >= 0.8]
  related_pairs_a <- related_pairs %>% group_by(sample_a) %>% summarise(r_a = paste(unique(sample_b), collapse=","))
  related_pairs_b <- related_pairs %>% group_by(sample_b) %>% summarise(r_b = paste(unique(sample_a), collapse=","))
  related_pairs <- merge(related_pairs_a, related_pairs_b, by.x="sample_a", by.y="sample_b", all=T) %>%
    mutate(related_samples=paste0(r_a,",",r_b)) %>% select(sample_a, related_samples) %>%
    separate_rows(related_samples, sep=",") %>% filter(related_samples != "NA") %>%
    group_by(sample_a) %>% summarise(related_samples = paste(unique(related_samples), collapse=",")) %>%
    rename(sample=sample_a)
  return(list(pairs=pairs_dt,related_samples=related_pairs,samples=samples_dt,ancestry=ancestry_dt))
}

read_alignment_QC <- function(project_id, data_path, somalier_path, canonical_chrs, minCov = 10,by_tranche=T,previous_rds=NA,force_update=F) {
  coverage <- list()
  samtools_stats <- list()
  readstats <- list()
  somalier <- list()
  summary_qc <- NULL
  sequencing_setup <- NULL
  old_rds <- NULL
  
  # If we need to update a previous RDS read it and clean previous applied filter columns
  if (!is.na(previous_rds)) {
    old_rds <- readRDS(previous_rds)
  }

  # When running by tranche, get a list of subfolders as tranche names
  tranches <- list()
  if (by_tranche) {
    tranche_name <- list.dirs(data_path, recursive = F,full.names = F)
    tranche_name <- tranche_name[ grepl("tranche", tranche_name) ]
    for (t in tranche_name) {
      tranche_path <- paste0(data_path,"/",t,"/qc_alignment")
      tranches[[t]] <- tranche_path
    }
  } else {
    tranches[[project_id]] <- paste0(data_path,"/qc_alignment")
  }

  # Load data from a previous RDS is present
  if (!is.null(old_rds)) {
    for (t in names(tranches)) {
      # if we want to force update then we need to remove new tranches that are aldready present
      if (force_update) {
        if (t %in% names(old_rds$coverage)) {
          old_rds$coverage[[t]] <- NULL
          old_rds$samtools_stats[[t]] <- NULL
          old_rds$readstats[[t]] <- NULL
          old_rds$somalier[[t]] <- NULL
          old_rds$summary_qc %>% filter(group != t)
          old_rds$sequencing_setup %>% filter(group != t)
        }
      } else {
        # If not forcing update than tranches already in the previous RDS will be skipped
        if (t %in% names(old_rds$coverage)) {
          tranches[[t]] <- NULL
        }
      }
    }
    # Store the relevant summary data so we can update the tables
    summary_qc <- old_rds$summary_qc %>% select(-QC_pass_status,-fail_reason)
    sequencing_setup <- old_rds$sequencing_setup
    coverage <- old_rds$coverage
    samtools_stats <- old_rds$samtools_stats
    readstats <- old_rds$readstats
    somalier <- old_rds$somalier
  } 

  # Cycle through tranches folders and read the QC data
  for (t in names(tranches)) {
    tranche_path <- tranches[[t]]
    coverage[[t]] <- read_mosdepth(paste0(tranche_path,"/coverage"), canonical_chrs,t,minCov)
    samtools_stats[[t]] <- list(
      idxstats = read_samtools_idxstats(paste0(tranche_path,"/mapflag_stats"), canonical_chrs, t),
      flagstats = read_samtools_flagstats(paste0(tranche_path,"/mapflag_stats"), t)
    )
    readstats[[t]] <- list(
      general = read_readstats(paste0(tranche_path,"/readstats"), t),
      s_per_run = read_samples_per_run(paste0(tranche_path,"/readstats"), t),
      r_per_sample = read_runs_per_sample(paste0(tranche_path,"/readstats"), t)
    )
    somalier[[t]] <- read_somalier(somalier_path, t)

    # Compile summary align QC table from single metrics tables
    message("### Create summary QC table ###")
    r_per_sample_wide <- as.data.table(readstats[[t]]$r_per_sample[, !"group"] %>% 
      pivot_wider(names_from = "metric", values_from = "n"))
    no_reads_contigs <- as.data.table(samtools_stats[[t]]$idxstats %>% 
      filter(n_mapped == 0) %>% group_by(sample) %>% 
      summarise(contigs_with_zero_reads = paste(contig, collapse=",")) %>%
      ungroup())
    list_of_dt <- list(coverage[[t]]$summarycov[, !"group"], 
                       samtools_stats[[t]]$flagstats[, .(sample, total, mapped, pct_mapped, primary, pct_primary, properly.paired, pct_properly.paired)], 
                       readstats[[t]]$general[, !"group"],
                       no_reads_contigs,
                       r_per_sample_wide, 
                       somalier[[t]]$ancestry[, .(sample,predicted_ancestry)], 
                       somalier[[t]]$samples[, .(sample,original_pedigree_sex,X_depth_mean,Y_depth_mean,scaled_X_depth,X_het)], 
                       somalier[[t]]$related_samples)
    merged_dt <- merge_tables(list_of_dt)
    merged_dt <- merged_dt %>% 
      mutate(inferred_sex = ifelse(
        (X_het < 5 & scaled_X_depth < 1.2) | (Y_depth_mean / X_depth_mean) >= 0.8 | (scaled_X_depth < 1.5 & Y_depth_mean / X_depth_mean > 0.5), "male", "female"))
    merged_dt$concordant_sex <- ifelse(merged_dt$inferred_sex == merged_dt$original_pedigree_sex, 1,0)
    merged_dt$group <- t
    summary_qc <- rbind(summary_qc, merged_dt)
    
    # Make a summary table for sequencing run setup
    sequencing_setup <- rbind(sequencing_setup,readstats[[t]]$s_per_run)
  }
  message("### ALIGN QC PROCESSING COMPLETE ###")
  return(
    list(
      coverage=coverage, 
      samtools_stats=samtools_stats,
      readstats=readstats,
      somalier=somalier,
      summary_qc=summary_qc,
      sequencing_setup=sequencing_setup
      )
    )
}

# optparse settings ----------------
option_list = list(
  make_option(c("-p", "--project_qc_dir"), action="store", default=NA, type='character',
              help="directory of the project QC data. Expects subfolders for tranches in it"),
  make_option(c("-s", "--somalier_tsv_dir"), action="store", default="./", type='character',
              help="directory where the TSV files from somalier ancestry/relate are stored"),
  make_option(c("-n", "--project_name"), action="store", default=NA, type='character',
              help="Project name. Used as subfolder to store data in the app"),
  make_option(c("-a", "--app_dir"), action="store", default=NA, type='character',
              help="Full path to the app directory"),
  make_option(c("-d", "--app_name"),  action="store", default=NA, type='character',
              help="name of the app to use for publishing to rsconnect"),
  make_option(c("-u", "--update_app"), action="store_true", default=FALSE,
              help="set this to update the shiny app"),
  make_option(c("-t", "--by_tranche"), action="store_true",
              help="set this when qc results are organized by tranche. Multiple sub-folders are then expected in project_qc_dir"),
  make_option(c("-f", "--align_fails"), action="store", default=NA, type='character',
              help="Text file containing filter to apply to align one per line. Samples selected by any of these are marked as FAIL"),
  make_option(c("-v", "--var_fails"), action="store", default=NA, type='character',
              help="Text file containing filter to apply to vars one per line. Samples selected by any of these are marked as FAIL"),
  make_option(c("-b", "--align_warnings"), action="store", default=NA, type='character',
              help="Text file containing filter to apply to align one per line. Samples selected by any of these are marked as WARNING"),
  make_option(c("-c", "--var_warnings"), action="store", default=NA, type='character',
              help="Text file containing filter to apply to vars one per line. Samples selected by any of these are marked as WARNING"),
  make_option(c("-x", "--force_update"), action="store_true", default=FALSE,
              help="Set this to force update data even if already present in the RDS"),
  make_option(c("-w", "--var_qc_rds"), action="store", default=NA, type='character',
              help="RDS file of previous variants qc data"),
  make_option(c("-z", "--align_qc_rds"), action="store", default=NA, type='character',
              help="RDS file of previous alignments qc data")                        
)
opt = parse_args(OptionParser(option_list=option_list))

# Check arguments -------------
project_dir = opt$project_qc_dir
project_name = opt$project_name
somalier_dir = opt$somalier_tsv_dir
app_dir = opt$app_dir
app_name = opt$app_name
align_fails = opt$align_fails
var_fails = opt$var_fails
align_warnings = opt$align_warnings
var_warnings = opt$var_warnings
by_tranche = ifelse(is.null(opt$by_tranche), FALSE, TRUE)
previous_var_rds = opt$var_qc_rds
previous_align_rds = opt$align_qc_rds
force_update = opt$force_update

if (opt$update_app) {
  if (is.na(app_dir) | is.na(app_name)) {
    message("ERROR - app dir and name must provided when update_app is active")
    quit(status = 1)
  }
}

if (!dir.exists(project_dir)) {
  message("ERROR - Directory ", project_dir, " does not exist")
  quit(status = 1)
}

# Read alignments QC and apply filters -------------
message("###############################")
message("###  ALIGN QC DATA PROCESS  ###")
message("###############################")
qc_align_data <- read_alignment_QC(project_name,project_dir, somalier_dir, canonical_chrs, min_cov, by_tranche, previous_var_rds, force_update)
qc_align_data$summary_qc$QC_pass_status <- "PASS"

message("### APPLY ALIGN FILTERS IF ANY ###")
qc_align_data$summary_qc <- evaluate_filters(qc_align_data$summary_qc, align_warnings, "WARNING")
qc_align_data$summary_qc <- evaluate_filters(qc_align_data$summary_qc, align_fails, "FAIL")

message("### STORE CURRENT TIMESTAMP ###")
qc_align_data$last_update <- timestamp

message("### SAVE ALIGN QC REPORTS ###")
saveRDS(qc_align_data, file = paste0(project_dir,"/qc_align.RDS"))
# check if qc_align_summary.tsv exists and if so, rename it with the timestamp
if (file.exists(paste0(project_dir,"/qc_align_summary.tsv"))) {
  file.rename(paste0(project_dir,"/qc_align_summary.tsv"), paste0(project_dir,"/qc_align_summary_",file_timestamp,".tsv"))
}
save_tsv(qc_align_data$summary_qc, paste0(project_dir,"/qc_align_summary.tsv"))

# Read variants QC and apply filters -------------
message("###############################")
message("###  VARS QC DATA PROCESS   ###")
message("###############################")
qc_vars_data <- read_vars_QC(project_name,project_dir, min_cov, min_qual, by_tranche, previous_var_rds, force_update)
if (!is.null(qc_vars_data$summary_qc)) {
  qc_vars_data$summary_qc$QC_pass_status <- "PASS"

  message("### APPLY VARS FILTERS IF ANY ###")
  qc_vars_data$summary_qc <- evaluate_filters(qc_vars_data$summary_qc, var_warnings, "WARNING")
  qc_vars_data$summary_qc <- evaluate_filters(qc_vars_data$summary_qc, var_fails, "FAIL")

  message("### STORE CURRENT TIMESTAMP ###")
  qc_vars_data$last_update <- timestamp

  message("### SAVE VARS QC REPORTS ###")
  saveRDS(qc_vars_data, file = paste0(project_dir,"/qc_vars.RDS"))

  # check if qc_vars_summary.tsv exists and if so, rename it with the timestamp
  if (file.exists(paste0(project_dir,"/qc_vars_summary.tsv"))) {
    file.rename(paste0(project_dir,"/qc_vars_summary.tsv"), paste0(project_dir,"/qc_vars_summary_",file_timestamp,".tsv"))
  }
  save_tsv(qc_vars_data$summary_qc, paste0(project_dir,"/qc_vars_summary.tsv"))
}

message("###############################")
message("###  MAKE PROJECT SUMMARY   ###")
message("###############################")
qc_align_fails <- qc_align_data$summary_qc %>% select(sample,group,FAIL_reason)
if (!is.null(qc_vars_data$summary_qc)) {
  qc_vars_fails <- qc_vars_data$summary_qc %>% select(sample,group,FAIL_reason)
} else {
   qc_vars_fails <- NULL
}
qc_details <- rbind(qc_align_fails, qc_vars_fails) 
qc_details <- qc_details %>% group_by(sample,group) %>% 
	summarise(FAIL_reason=paste(unique(FAIL_reason), collapse=";")) %>%
	mutate(FAIL_reason=gsub("(;PASS|PASS;|;PASS;)","",FAIL_reason)) %>%
  ungroup()
qc_details <- qc_details %>% group_by(group,FAIL_reason) %>% 
	summarise(N=n()) %>% ungroup()

# in each group in qc_details, if FAIL_reason do not contain "PASS" then add a row with PASS and N=0
for (g in unique(qc_details$group)) {
  if (!"PASS" %in% qc_details$FAIL_reason[qc_details$group == g]) {
    qc_details <- rbind(qc_details, data.frame(group=g,FAIL_reason="PASS",N=0))
  }
}

qc_details <- qc_details %>% pivot_wider(names_from="FAIL_reason", values_from="N", values_fill=0)
qc_details <- qc_details %>% relocate(PASS, .after=group)
if (ncol(qc_details) == 2) {
  qc_details$FAIL <- 0
} else {
  qc_details$FAIL <- rowSums(qc_details[,3:ncol(qc_details)])
}
qc_details <- qc_details %>% relocate(FAIL, .after=PASS)
qc_details$tot_samples <- qc_details$PASS + qc_details$FAIL
qc_details <- qc_details %>% relocate(tot_samples, .after=group)

header <- c("# FAIL REASONS")
if (ncol(qc_details) > 4) {
  x <- 0
  for (n in 5:ncol(qc_details)) {
    x <- x + 1
    filter_tag <- colnames(qc_details)[n]
    id <- paste0("F",x)
    colnames(qc_details)[n] <- id
    header <- c(header, paste0("# ", id, ": ", filter_tag))
  }
}
qc_details <- rbind(as.data.frame(qc_details), c("Total", colSums(qc_details[,2:ncol(qc_details)])))
qc_details <- qc_details %>% mutate(across(2:ncol(.), as.numeric))
qc_details$fraction_FAIL <- round(qc_details$FAIL / qc_details$tot_samples,3)
qc_details <- qc_details %>% relocate(fraction_FAIL, .after=FAIL)
write(header, file=paste0(project_dir,"/project_qc_summary.tsv"), append=F)
write.table(qc_details, file=paste0(project_dir,"/project_qc_summary.tsv"), sep="\t", row.names=F, quote=F, append=T)

# Update shiny app data if requested ------------
if (opt$update_app) {
  message("#########################")
  message("###  UPDATE APP DATA  ###")
  message("#########################")
  app_project_dir = paste0(app_dir,"/data/",project_name)
  dir.create(app_project_dir,recursive = T)
  file.copy(paste0(project_dir,"/qc_align.RDS"),paste0(app_project_dir, "/qc_align.RDS"),overwrite = T)
  file.copy(paste0(project_dir,"/qc_vars.RDS"),paste0(app_project_dir, "/qc_vars.RDS"),overwrite = T)
  file.copy(paste0(project_dir,"/project_qc_summary.tsv"),paste0(app_project_dir, "/project_qc_summary.tsv"),overwrite = T)

  server = Sys.getenv("RSCONNECT_SERVER")
  key = Sys.getenv("RSCONNECT_APIKEY")
  addServer(server, name="connect-server")
  connectApiUser(account = "deploy_account", apiKey = key, quiet = FALSE, server="connect-server")
  deployApp(appDir=app_dir, appName=app_name, account="deploy_account", server="connect-server", logLevel="normal",upload=T, forceUpdate = TRUE)
}

message("### All done! ###")