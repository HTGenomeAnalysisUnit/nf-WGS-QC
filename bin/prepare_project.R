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
canonical_chrs <- paste("chr",c(1:22,"X","Y","M","MT"), sep="")

# Helper functions -----------------
`%nin%` <- Negate('%in%')

save_tsv <- function(df, outfile) {
  write.table(df, file = outfile, sep="\t", quote=F, row.names = F)
}

evaluate_filters <- function(df, filters_file) {
  if (!is.na(filters_file)) {
    df$QC_pass_status <- "PASS"
    df$fail_reason <- "PASS"
    for (filter_string in scan(filters_file, what="", sep="\n")) {
      df <- df %>% 
        mutate(
          QC_pass_status = ifelse(
            eval_tidy(parse_expr(filter_string)), 
            "FAIL", 
            QC_pass_status
          ),
          fail_reason = ifelse(
            eval_tidy(parse_expr(filter_string)), 
            paste(fail_reason,filter_string,sep=";"), 
            fail_reason
          )
        ) 
    }
    df <- df %>% mutate(
        fail_reason=ifelse(
          fail_reason=="PASS",
          "PASS",
          gsub("PASS;","",fail_reason)
        )
      )  
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
  res <- list(by_sample=list(), all_samples=NULL)
  
  files <- list.files(path, pattern = "varsPerContig.tsv")
  for (f in files) {
    sample_name <- gsub("\\.varsPerContig\\.tsv","", f)
    dt <- fread(paste0(path, "/", f), header=F, col.names=c("contig","size","n_vars"))
    idx <- which(dt$contig %in% contigs)
    dt <- dt[idx]
    dt$vars_per_Mb <- dt$n_vars / (dt$size / 1000000)
    dt$sample <- sample_name
    res$by_sample[[sample_name]] <- dt
    res$all_samples <- rbind(res$all_samples, dt)
  }
  return(res)
}

read_stats_vars <- function(path, minDP, minQUAL) {
  res <- list(by_sample=list(), DP=NULL, QUAL=NULL, summary=NULL)
  
  files <- list.files(path, pattern = "vars.stats")
  for (f in files) {
    sample_name <- gsub("\\.vars\\.stats","", f)
    sdata <- readVCFstat(paste0(path, "/", f), fileID = sample_name)
    res$by_sample[[sample_name]] <- sdata
    
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
  files <- list.files(path, pattern = "cont.tsv")
  for (f in files) {
    dt <- fread(paste0(path, "/", f), header=T, check.names=T)
    setnames(dt, "X.SAMPLE", "sample")
    charr_stats <- list(mean=mean(dt$CHARR), sd=sd(dt$CHARR))
    het_rate_stats <- list(mean=mean(dt$HETEROZYGOSITY_RATE), sd=sd(dt$HETEROZYGOSITY_RATE))
    dt$charr_sd_from_mean <- round(abs(dt$CHARR - charr_stats$mean) / charr_stats$sd, 3)
    dt$hetrate_sd_from_mean <- round(abs(dt$CHARR - het_rate_stats$mean) / het_rate_stats$sd, 3)
  }
  return(dt)
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
    stats[[t]] <- read_stats_vars(var_stats_path, minDP, minQUAL)
    vars_per_contig[[t]] <- read_contig_vars(var_stats_path, canonical_chrs)
    contamination[[t]] <- read_contamination_data(contamination_data_path)

    list_of_dt <- list(
      stats[[t]]$summary,
      contamination[[t]][,.(sample, CHARR, HETEROZYGOSITY_RATE,charr_sd_from_mean,hetrate_sd_from_mean)]
    )
    merged_dt <- merge_tables(list_of_dt)
    merged_dt$group <- t
    summary_qc <- rbind(summary_qc, merged_dt)
  }
  
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
  files <- list.files(path, pattern="global.dist.txt")
  for (f in files) {
    message("Processing ", f)
    sample_name <- gsub("\\.mosdepth\\.global\\.dist\\.txt","",f)
    dt <- fread(cmd = paste0("grep total ", path, "/", f), header=F, 
                col.names=c("min_cov","value"), drop=c(1))
    dt[, sample := sample_name ]
    cumcov_dt <- rbind(cumcov_dt, dt)
    mc <- dt[min_cov == minCov][, !"min_cov"]
    setnames(mc, "value", paste0("pct_bases_min_",minCov))
    mincov_dt <- rbind(mincov_dt, mc)
  }
  cumcov_dt$group <- group_name
  
  meancov_dt <- NULL
  summary_path <- paste0(path, "/summary")
  files <- list.files(summary_path, pattern="*.summary.txt")
  for (f in files) {
    message("Processing ", f)
    sample_name <- gsub("\\.mosdepth\\.summary\\.txt","",f)
    dt <- fread(paste0(summary_path, "/", f), header=T, 
                select=c("chrom","mean","min","max"))
    setnames(dt, old = c("mean","min","max"), new = c("mean_cov","min_cov","max_cov"))
    contigs_idx <- c(which(dt$chrom %in% canonical_chrs), which(dt$chrom == "total"))
    dt <- dt[contigs_idx]
    dt <- dt[, sample := sample_name ]
    meancov_dt <- rbind(meancov_dt, dt)
  }
  meancov_dt$group <- group_name
  summary_cov <- merge(meancov_dt[chrom == "total"], mincov_dt, by="sample", all.x=T)
  return(list(cumcov=cumcov_dt, meancov=meancov_dt,summarycov=summary_cov))
}

read_samtools_idxstats <- function(path, canonical_chrs, group_name, decoy_string="decoy") {
  res_dt <- NULL
  files <- list.files(path, pattern="\\.idxstat")
  for (f in files) {
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
  files <- list.files(path, pattern="\\.flagstat")
  for (f in files) {
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
  files <- list.files(path, pattern="\\.per_sample_run\\.tsv")
  for (f in files) {
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
  pairs_dt <- fread(paste0(path, "/",group_name, ".pairs.tsv"), header=T)
  setnames(pairs_dt, "#sample_a", "sample_a")
  pairs_dt$group <- group_name
  samples_dt <- fread(paste0(path, "/",group_name, ".samples.tsv"), header=T, 
                      select=c("sample_id","original_pedigree_sex","X_depth_mean","Y_depth_mean","n_hom_ref","n_het","n_hom_alt"))
  samples_dt$group <- group_name
  setnames(samples_dt, "sample_id", "sample")
  ancestry_dt <- fread(paste0(group_name, ".somalier-ancestry.tsv"), header=T)
  ancestry_dt$group <- group_name
  ancestry_dt$given_ancestry[ancestry_dt$given_ancestry == ""] <- "External_sample"
  setnames(ancestry_dt, "#sample_id", "sample")
  
  related_pairs <- pairs_dt[relatedness >= 0.5]
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

    # Compile summary aling QC table from single metrics tables
    r_per_sample_wide <- as.data.table(readstats[[t]]$r_per_sample[, !"group"] %>% 
      pivot_wider(names_from = "metric", values_from = "n"))
    list_of_dt <- list(coverage[[t]]$summarycov[, !"group"], 
                       samtools_stats[[t]]$flagstats[, .(sample, total, mapped, pct_mapped, primary, pct_primary, properly.paired, pct_properly.paired)], 
                       readstats[[t]]$general[, !"group"], 
                       r_per_sample_wide, 
                       somalier[[t]]$ancestry[, .(sample,predicted_ancestry)], 
                       somalier[[t]]$samples[, .(sample,original_pedigree_sex,X_depth_mean,Y_depth_mean)], 
                       somalier[[t]]$related_samples)
    merged_dt <- merge_tables(list_of_dt)
    merged_dt$inferred_sex <- ifelse(
      (merged_dt$Y_depth_mean / merged_dt$X_depth_mean) > 0.8 & (merged_dt$Y_depth_mean / merged_dt$X_depth_mean) < 1.2,
      "male", "female")
    merged_dt$concordant_sex <- ifelse(merged_dt$inferred_sex == merged_dt$original_pedigree_sex, 1,0)
    merged_dt$group <- t
    summary_qc <- rbind(summary_qc, merged_dt)
    
    # Make a summary table for sequencing run setup
    sequencing_setup <- rbind(sequencing_setup,readstats[[t]]$s_per_run)
  }

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
  make_option(c("-f", "--align_filters"), action="store", default=NA, type='character',
              help="Text file containing filter to apply to align one per line. Samples selected by any of these are marked as FAIL"),
  make_option(c("-v", "--var_filters"), action="store", default=NA, type='character',
              help="Text file containing filter to apply to vars one per line. Samples selected by any of these are marked as FAIL"),
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
align_filters = opt$align_filters
var_filters = opt$var_filters
by_tranche = opt$by_tranche
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

# Read alignements QC and apply filters -------------
message("Reading align QC data...")
qc_align_data <- read_alignment_QC(p,project_dir, somalier_dir, canonical_chrs, min_cov, by_tranche, previous_var_rds, force_update)
qc_align_data$summary_qc$QC_pass_status <- "NOT_APPLIED"

qc_align_data$summary_qc <- evaluate_filters(qc_align_data$summary_qc, align_filters)

saveRDS(qc_align_data, file = paste0(project_dir,"/qc_align.RDS"))
save_tsv(qc_align_data$summary_qc, paste0(project_dir,"/qc_align_summary.tsv"))

# Read variants QC and apply filters -------------
message("Reading variants QC data...")
qc_vars_data <- read_vars_QC(p,project_dir, min_cov, min_qual, by_tranche, previous_var_rds, force_update)
qc_vars_data$summary_qc$QC_pass_status <- "NOT_APPLIED"

qc_vars_data$summary_qc <- evaluate_filters(qc_vars_data$summary_qc, var_filters)

saveRDS(qc_vars_data, file = paste0(project_dir,"/qc_vars.RDS"))
save_tsv(qc_vars_data$summary_qc, paste0(project_dir,"/qc_vars_summary.tsv"))

qc_align_merge <- qc_align_data$summary_qc %>% select(sample,group,QC_pass_status)
qc_vars_merge <- qc_vars_data$summary_qc %>% select(sample,QC_pass_status)
setnames(qc_align_merge, "QC_pass_status", "align_pass_status")
setnames(qc_vars_merge, "QC_pass_status", "vars_pass_status")
qc_merged <- merge(qc_align_merge, qc_vars_merge, by="sample", all.x=T)
qc_merged %>% mutate(pass_status=ifelse(align_pass_status == "FAIL" | vars_pass_status == "FAIL", "FAIL", "PASS"))

qc_merged_counts <- qc_merged %>% group_by(group) %>%
  summarise(n_samples=n(), n_pass=sum(pass_status == "PASS"), n_fail=sum(pass_status == "FAIL")) %>%
  mutate(fraction_fail=n_fail/n_samples)
qc_merged_counts <- c("total", 
  sum(qc_merged_counts$n_samples),
  sum(qc_merged_counts$n_pass),
  sum(qc_merged_counts$n_fail),
  sum(qc_merged_counts$n_fail)/sum(qc_merged_counts$n_samples))
save_tsv(qc_merged_counts, paste0(project_dir,"/project_qc_summary.tsv"))

# Update shiny app data if requested ------------
if (opt$update_app) {
  message("Update app data")
  app_project_dir = paste0(app_dir,"/data/",project_name)
  dir.create(app_project_dir,recursive = T)
  file.copy(paste0(project_dir,"/qc_align.RDS"),paste0(app_project_dir, "/qc_align.RDS"),overwrite = T)
  file.copy(paste0(project_dir,"/qc_vars.RDS"),paste0(app_project_dir, "/qc_vars.RDS"),overwrite = T)
  file.copy(paste0(project_dir,"/project_qc_summary.tsv"),paste0(app_project_dir, "/project_qc_summary.tsv"),overwrite = T)

  server = Sys.getenv("RSCONNECT_SERVER")
  key = Sys.getenv("RSCONNECT_APIKEY")
  connectApiUser(account = "deploy_account", apiKey = key, quiet = FALSE, server=server)
  deployApp(appName=app_name, account="deploy_account", logLevel="normal",upload=T)
}

message("All done!")