library(data.table)
library(ggplot2)
library(ggrepel)
set_theme(theme_classic())
set.seed(150799)
output_dir="PEPITOPE/s7_screen_results/answering_JV_observations/missing_barcodes/"
if(!dir.exists(output_dir)){dir.create(output_dir)}

# Load data
counts_screen_1=fread("PEPITOPE/s7_screen_results/output_screen_1/manual_qc_rows.txt")
counts_screen_2a=fread("PEPITOPE/s7_screen_results/output_screen_2a/manual_qc_rows.txt")

# # Grab TAA barcodes
taa_big=fread("PEPITOPE/s6_library_runs/shared_by_JV/large_TAA.csv")
setnames(taa_big, "BARCODE", "barcode")
taa_big=taa_big[!is.na(LENGTH)]
setnames(taa_big, "GENE NAME 2", "gene_name")
setnames(taa_big, "MINIGENE SEQ", "tiled")
taa_big[, pep_type:= "TAA"]
taa_big[, mut_id:= "TAA"]
taa_big[, pep_id:= "TAA"]
taa_big[, sample_id:= "TAA_big"]
taa_big[, patient_id:= "TAA_big"]

# #  Pat neoantigens
all_constructs=fread("PEPITOPE/s5_create_libraries/new_rna_library_all_patients_duplo.txt")
setnames(all_constructs, "MINIGENE", "tiled")
setnames(all_constructs, "BARCODE", "barcode")
pat_dt=all_constructs[, .(barcode, replicate, patient_id, gene_name, mut_id, pep_id, pep_type)]

# # create patients + TAA big
tot_barcode_anno_big=rbind(
  taa_big[, .(barcode, replicate="1", patient_id, gene_name, mut_id, pep_id, pep_type)], 
  pat_dt)

tot_barcode_anno_big[, gene_id:= ifelse(pep_id=="TAA", gene_name,
                                        ifelse(gene_name=="", gsub("\\.fs", "",mut_id),gene_name))]

tot_barcode_anno_big[, variable_id:=paste(gene_id,pep_type,replicate,sep = "_")]
fwrite(tot_barcode_anno_big, paste0(output_dir, "tot_barcode.txt"))


check_incomplete_cases <- function(counts_screen, pat,
                                   goal  = tot_barcode_anno_big,
                                   grain = "minigene") {
  
  ## 1. GOAL: expected constructs for this patient (+ shared TAA library)
  goal_pat <- goal[patient_id %in% c("TAA_big", pat)]
  # construct unit:
  #   minigene -> pep_id for neoantigens (duplo = 2 barcodes); gene_name for TAA (pep_id == "TAA")
  #   gene     -> gene_name (pools all mutations of a gene: TTN alt = 94 barcodes as one unit)
  goal_pat[, unit := if (grain == "minigene")
    fifelse(pep_type == "TAA", gene_name, pep_id)
    else gene_name]
  goal_dt <- goal_pat[, .(n_expected = .N), by = .(unit, mut_profile = pep_type)]
  
  ## 2. ACTUAL: counts -> long, annotated by BARCODE (guide == barcode). No gene-string parsing.
  id_cols     <- c("guide", "gene", "guide_type")
  sample_cols <- setdiff(names(counts_screen), id_cols)
  mcounts <- melt(counts_screen, id.vars = id_cols, measure.vars = sample_cols,
                  variable.name = "variable", value.name = "value")
  mcounts <- merge(mcounts,
                   goal_pat[, .(guide = barcode, unit, mut_profile = pep_type)],
                   by = "guide", all.y=T)          # inner join -> keeps only this patient's + TAA barcodes
  
  ## optional QC: per-barcode count distribution
  ggplot(mcounts, aes(variable, value + 1)) +
    geom_boxplot(aes(fill = mut_profile), outlier.size = .1) +
    geom_hline(yintercept = c(100, 1000), lty = "dashed", lwd = .3) +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    scale_fill_manual(values = c(alt = "red3", ref = "grey60", TAA = "gold"))
  ggsave(paste0(output_dir, pat, "_barcode_counts.pdf"), width = 6, height = 4)
  
  ## 3. detected barcodes per unit x sample
  det <- mcounts[, .(n_detected = sum(value > 0)), by = .(unit, mut_profile, variable)]
  
  ## 4. full grid (unit x sample) so units absent everywhere stay in the denominator
  samples <- unique(mcounts$variable)
  full <- goal_dt[, .(variable = samples), by = .(unit, mut_profile, n_expected)]
  full <- merge(full, det, by = c("unit", "mut_profile", "variable"), all.x = TRUE)
  full[is.na(n_detected), n_detected := 0L]
  full[, status := fifelse(n_detected == 0L,        "dropout",   # total loss (0 barcodes)
                           fifelse(n_detected <  n_expected, "partial",   # some barcodes missing
                                   "complete"))]
  fwrite(full, paste0(output_dir, pat, "_completeness_per_unit.txt"))
  
  ## 5. per-sample summary
  comp <- full[, .(tot      = .N,
                   complete = sum(status == "complete"),
                   partial  = sum(status == "partial"),
                   dropout  = sum(status == "dropout")),
               by = .(variable, mut_profile)]
  comp[, pct_complete := 100 * complete / tot]
  fwrite(comp, paste0(output_dir, pat, "_completeness_stats.txt"))
  
  p1 <- ggplot(comp[variable != "unmatched"], aes(variable, pct_complete)) +
    geom_point(aes(col = mut_profile), position = position_dodge(width=.3)) + ylab("% complete (all barcodes present)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    scale_color_manual(values = c(alt = "red3", ref = "black", TAA = "gold"))
  
  p2 <- ggplot(comp[variable != "unmatched"], aes(variable, tot - complete)) +
    geom_point(aes(col = mut_profile), position = position_dodge(width=.3)) + ylab("n incomplete units") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    scale_color_manual(values = c(alt = "red3", ref = "black", TAA = "gold"))
  
  cowplot::plot_grid(p1, p2)
  ggsave(paste0(output_dir, pat, "_viz_incomplete_cases.pdf"), width = 9, height = 4)
  
  # IF TAKING 100 as a measure
  
  ## 3. detected barcodes per unit x sample
  det <- mcounts[, .(n_detected = sum(value > 100)), by = .(unit, mut_profile, variable)]
  
  ## 4. full grid (unit x sample) so units absent everywhere stay in the denominator
  samples <- unique(mcounts$variable)
  full <- goal_dt[, .(variable = samples), by = .(unit, mut_profile, n_expected)]
  full <- merge(full, det, by = c("unit", "mut_profile", "variable"), all.x = TRUE)
  full[is.na(n_detected), n_detected := 0L]
  full[, status := fifelse(n_detected == 0L,        "dropout",   # total loss (0 barcodes)
                           fifelse(n_detected <  n_expected, "partial",   # some barcodes missing
                                   "complete"))]
  fwrite(full, paste0(output_dir, pat, "_completeness_per_unit.txt"))
  
  ## 5. per-sample summary
  comp <- full[, .(tot      = .N,
                   complete = sum(status == "complete"),
                   partial  = sum(status == "partial"),
                   dropout  = sum(status == "dropout")),
               by = .(variable, mut_profile)]
  comp[, pct_complete := 100 * complete / tot]
  fwrite(comp, paste0(output_dir, pat, "_completeness_stats_taking_100_reads_min.txt"))
  
  p1 <- ggplot(comp[variable != "unmatched"], aes(variable, pct_complete)) +
    geom_point(aes(col = mut_profile), position = position_dodge(width=.3)) + ylab("% complete (all barcodes present)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    scale_color_manual(values = c(alt = "red3", ref = "black", TAA = "gold"))
  
  p2 <- ggplot(comp[variable != "unmatched"], aes(variable, tot - complete)) +
    geom_point(aes(col = mut_profile), position = position_dodge(width=.3)) + ylab("n incomplete units") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    scale_color_manual(values = c(alt = "red3", ref = "black", TAA = "gold"))
  
  cowplot::plot_grid(p1, p2)
  ggsave(paste0(output_dir, pat, "_viz_incomplete_cases_taking_100_reads_min.pdf"), width = 9, height = 4)
  
  invisible(full)
}

check_incomplete_cases(counts_screen_1,  "P1")
check_incomplete_cases(counts_screen_2a, "P5")
