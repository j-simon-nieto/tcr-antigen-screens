library(pepitope)
library(data.table)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
library(readxl)
library(Biostrings)
library(DESeq2)
set.seed(150799)
theme_set(theme_classic())

output_dir="PEPITOPE/s7_NEW_screen_results/output_screen_2a/updated_screen_results/"
# updates to the logfc thresholds in deseq2
# run on 18/Sept/26: fix fusion proteins
if(!dir.exists(output_dir)){
  dir.create(output_dir)}
source("PEPITOPE/s7_NEW_screen_results/palette_screens.R")

# Grabbing TAA big library for barcodes

# 0/ Load files
# ================================================

# # results from QC
stats_dt=fread(paste0(output_dir, "../manual_qc_stats.txt"))
counts_dt=fread(paste0(output_dir, "../manual_qc_counts.txt"))
rows_dt=fread(paste0(output_dir, "../manual_qc_rows.txt"))

# 1/ Re-format and adapt DESEQ2
# ================================================

# create central counts and design
counts_mat= as.matrix(as.data.frame(counts_dt)[,3:ncol(counts_dt)])
rownames(counts_mat)=counts_dt$guide

design_dt=data.table(sample_id=colnames(counts_mat))
design_dt[, rep:= stringr::str_split_i(sample_id, " ", -1)]
design_dt[, origin:= stringr::str_split_i(sample_id, " ", 1)]

# --------------- RENAMING STARTS ----------------------
# ------------------------------------------------------
# # revisiting on 8th Oct to include the reference pairs to some TAAs:
hotspots=fread("PEPITOPE/s7_NEW_screen_results/screen_interpretation/Large_TAA2.0_hotspots.csv")
hotspots=hotspots[,.(gene_id=pep_id, mutation_profile= "TAA", type_TAA=type, ref_key)]
hotspots[, ref_key:= sub(".$", "", ref_key)]
# # model antigens: mut or ref appended as indication
model_antigens=data.table(gene_id=scan("PEPITOPE/s7_NEW_screen_results/screen_interpretation/Model_antigens.txt", 
                                       what="characters"),mutation_profile= "TAA")
model_antigens=model_antigens[grepl("mut|wt", gene_id),]
model_antigens[, type_TAA:= ifelse(grepl("mut$", gene_id), "alt", "ref")]
model_antigens[, ref_key:= sub("mut$|wt$", "", gene_id)]

taa_wt_ref_dt=rbindlist(list(hotspots, model_antigens))
# loading new info ends here

setnames(rows_dt, "guide", "barcode")
rows_dt <- rows_dt[match(rownames(counts_mat), rows_dt$barcode), ]

# # basic columns
rows_dt[, gene_replicate:=stringr::str_split_i(gene, "_", -1)]
rows_dt[, mutation_profile:=stringr::str_split_i(gene, "_", -2)]
rows_dt[, construct := sub("_[^_]+_[^_]+$", "", gene)] 

rows_dt[, gene_name := fifelse(mutation_profile == "TAA",
                               sub("_TAA$", "", construct),   
                               sub("_.*$", "", construct))]  
rows_dt[, pep_id := fifelse(mutation_profile == "TAA", NA_character_,
                            sub("^[^_]*_", "", construct))]
rows_dt[, is_fusion := grepl("--", pep_id, fixed = TRUE)]
rows_dt[, gene_id := fcase(mutation_profile == "TAA", gene_name,
                           default = pep_id)]

.locus <- function(x) sub("(?<=[0-9])(fs|[A-Za-z]+)[*]?$", "", sub("-[0-9]+$", "", x), perl = TRUE)

# alts declare their key(s)
rows_dt[, ref_key_3p := fcase(
  mutation_profile != "alt", NA_character_,
  is_fusion, sub("\\.fs$", "", sub("-[0-9]+$", "", sub("^.*--", "", pep_id))),
  default = .locus(pep_id))]
rows_dt[, ref_key_5p := fifelse(mutation_profile == "alt" & is_fusion,
                                sub("--.*$", "", pep_id), NA_character_)]

alt_keys <- unique(stats::na.omit(c(rows_dt$ref_key_3p, rows_dt$ref_key_5p)))

# refs adopt whichever of their own candidate keys an alt is asking for
rows_dt[mutation_profile == "ref", ref_key := {
  cand <- list(pep_id, .locus(pep_id), gene_name, .locus(gene_name))
  out  <- rep(NA_character_, .N)
  for (cc in cand) out <- fifelse(is.na(out) & cc %in% alt_keys, cc, out)
  fcoalesce(out, .locus(pep_id))
}]
rows_dt[mutation_profile == "alt", ref_key := ref_key_3p]  

# # 8th Oct: add new data on TAA alt/wt
reclassified_taa_ids=taa_wt_ref_dt$gene_id
rows_dt[taa_wt_ref_dt[gene_id %in% reclassified_taa_ids], on = "gene_id",
        `:=`(type_TAA = i.type_TAA, ref_key = i.ref_key)]

rows_dt[, mutation_profile := fcase(
  mutation_profile %in% c("alt","ref"),paste0("NeoAg_", mutation_profile),
  mutation_profile == "TAA" & type_TAA %in% "ref","TAA_ref",
  mutation_profile == "TAA" & grepl("CMV|EBV|HPV|FLU|SARS", gene_id), "ViralAg",
  mutation_profile == "TAA","TAA",
  default = NA_character_)]

# --------------- end of renaming efforts ----------------------
# --------------------------------------------------------------

design_dt=design_dt[sample_id!="unmatched",]
counts_mat=counts_mat[, colnames(counts_mat) %in% design_dt$sample_id]

stopifnot(identical(rows_dt$barcode, rownames(counts_mat)) & identical(design_dt$sample_id, colnames(counts_mat)))

dset <- DESeqDataSetFromMatrix(counts_mat, design_dt, ~ origin, rowData=rows_dt)

# Deseq function
screen_calc = function(dset, comparisons, min_count=30) {
  if (length(unique(dset$patient)) > 1)
    stop("The 'patient' column can not span more than one value")
  if (is.numeric(dset$rep))
    dset$rep = factor(dset$rep)
  
  # eset = DESeq2::DESeqDataSet(dset, ~ rep + origin) # removing because the replicates are not paired!!!
  eset = DESeq2::DESeqDataSet(dset, ~ origin)
  
  # pre-filter: drop barcodes that never reach min_count in ANY sample.
  # Applied before size factors so normalization ignores dropped-out/contaminant
  # barcodes; global (max over all samples) so padj stays comparable across contrasts.
  keep = matrixStats::rowMaxs(as.matrix(SummarizedExperiment::assay(eset))) >= min_count
  message(sprintf(
    "screen_calc: keeping %d/%d barcodes with >= %d counts in at least one sample (dropped %d)",
    sum(keep), length(keep), min_count, sum(!keep)))
  eset = eset[keep, ]
  
  eset$origin = factor(make.names(eset$origin))
  eset$rep = factor(eset$rep)
  DESeq2::sizeFactors(eset) = colSums(assay(eset)) / max(colSums(assay(eset)))
  mod = DESeq2::DESeq(eset, fitType="local")
  
  get_result = function(comp) {
    DESeq2::results(mod, contrast=c("origin", comp)) |>
      # lfcThreshold = 0.1, # significance tested as stronger than 0.25
      # alpha = .1, # this is default
      # altHypothesis = "greaterAbs"
      as.data.frame() |>
      tibble::rownames_to_column("barcode") |>
      as_tibble() |>
      filter(!is.na(log2FoldChange)) |>
      arrange(padj, pvalue) |>
      left_join(as.data.frame(rowData(eset)), by=join_by(barcode))
  }
  
  if (is.list(comparisons)) {
    lapply(comparisons, get_result) |>
      setNames(sapply(comparisons, paste, collapse=" vs "))
  } else {
    get_result(comparisons)
  }}

# Run DESEQ

results=screen_calc(dset, comparisons=list(c("UT", "B"),
                                           c("Single", "UT"),
                                           c("Cluster", "UT"),
                                           c("Single", "B"),
                                           c("Cluster", "B")))

results_agg=rbindlist(lapply(results, setDT), idcol = "comparison")
fwrite(results_agg, paste0(output_dir, "screen2a_results.txt"))

ggplot(results_agg[baseMean>30], aes(log2FoldChange, -log10(padj)))+geom_point(aes(col=mutation_profile))+
  ggrepel::geom_text_repel(data=results_agg[padj<1e-10],aes(label=gene_name))+
  facet_wrap(~comparison)+
  scale_color_manual(values=palette_mutprof)
ggsave(paste0(output_dir, "volcano_plot_overview.pdf"), width = 12, height=8)

ggplot(results_agg[baseMean>30], aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_agg[padj<1e-5],aes(label=gene_name), size=2.2)+
  facet_wrap(~comparison)+scale_x_log10()+
  geom_point(shape=1, data=results_agg[padj<0.05], size=1)+
  scale_color_manual(values=palette_mutprof)
ggsave(paste0(output_dir, "MA_like_plot_overview.pdf"), width = 12, height=8)

ggplot(results_agg[baseMean>100 & (comparison%in% c("Cluster vs UT", "Single vs UT"))], aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_agg[baseMean>100 &padj<0.05& comparison%in% c("Cluster vs UT", "Single vs UT")],
                           aes(label=gene_id), size=2.5, max.overlaps = 15)+
  facet_wrap(~comparison)+scale_x_log10()+
  geom_point(shape=1, data=results_agg[baseMean>100 &padj<0.05& comparison%in% c("Cluster vs UT", "Single vs UT")], size=1)+
  scale_color_manual(values=palette_mutprof)
ggsave(paste0(output_dir, "MA_like_plot_clusters_and_singlets_vs_UT.pdf"), width = 11, height=4)

ggplot(results_agg[baseMean>100 & grepl("vs UT", comparison)], aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_agg[baseMean>100 &padj<1e-5 & grepl("vs UT", comparison)],
                           aes(label=gene_id), size=2.5)+
  facet_wrap(~comparison)+scale_x_log10()+
  geom_point(shape=1, data=results_agg[baseMean>100 &padj<0.05 & grepl("vs UT", comparison)], size=1)+
  scale_color_manual(values=palette_mutprof)+
  geom_hline(yintercept = 0)+
  theme(legend.position = "bottom", axis.title = element_text(size=14))
ggsave(paste0(output_dir, "Clus_vs_UT_and_Singlet_vs_UT.pdf"), width = 7, height=4)

comp_db_sg_ut=merge(results$`Cluster vs UT`, results$`Single vs UT`, by=names(rows_dt)[!grepl(" ", names(rows_dt))], suffixes=c("_db", "_sg"))

ggplot(comp_db_sg_ut[baseMean_db>30 | baseMean_sg>30], aes(stat_sg, stat_db))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=comp_db_sg_ut[(baseMean_db>30 | baseMean_sg>30) & 
                                                (padj_db<0.05|padj_sg<0.05)],
                           aes(label=gene_id), size=3, max.overlaps = 10)+
  scale_color_manual(values=palette_mutprof)+
  labs(y="Stat: Cluster TCRs vs UT", x="Stat: Singlet TCRs vs UT")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")
ggsave(paste0(output_dir, "Clus_vs_UT__Singlet_vs_UT.pdf"), width = 5, height=5)

# # Analysis of top hits:
expand_for_plot <- function(dt) {
  alts <- copy(dt[mutation_profile != "ref"])
  alts[, plot_id := gene_id]
  
  refs <- copy(dt[mutation_profile == "ref"])
  map  <- unique(alts[!is.na(ref_key), .(ref_key, plot_id)])
  refs <- refs[map, on = "ref_key", allow.cartesian = TRUE]
  
  rbind(alts, refs, use.names = TRUE, fill = TRUE)
}

plot_top_dropouts = function(res_dt, top_n = 40, filename_pdf) {
  
  res_dt <- res_dt[baseMean > 30]
  
  # thresholds from the UNEXPANDED table -- refs are duplicated after expansion
  # and would otherwise be double-counted in mean/sd
  sd_x <- sd(res_dt$stat); mean_x <- mean(res_dt$stat)
  
  pd <- expand_for_plot(res_dt)
  pd[, min_stat := min(stat[mutation_profile != "ref"]), by = plot_id]  # rank on alts only
  pd=pd[!is.na(min_stat)]
  
  setorder(pd, min_stat)
  top <- unique(pd$plot_id)[seq_len(min(top_n, uniqueN(pd$plot_id)))]
  
  ggplot(pd[plot_id %in% top],
         aes(reorder(plot_id, -min_stat), stat, col = mutation_profile)) +
    geom_point(size = 0.8) + coord_flip() +
    ggbeeswarm::geom_beeswarm(cex = .15) +
    geom_hline(yintercept = c(mean_x - sd_x, mean_x, mean_x + sd_x),
               lty = c("dashed", "solid", "dashed"), lwd = c(.3, .5, .3)) +
    scale_color_manual(values = palette_mutprof) +
    geom_point(shape = 1, data = pd[plot_id %in% top & padj < 0.05],
               size = 3, col = "black") +
    theme(legend.position = "none", axis.title = element_text(size = 14)) +
    scale_y_reverse() + xlab("") + ylab("Drop-out confidence (T-stat)")
  ggsave(paste0(output_dir, filename_pdf), width = 6, height = 6)
}

res_clus=results$`Cluster vs UT`
plot_top_dropouts(res_clus, filename_pdf="top40_cluster_tcr_vs_ut_dropouts.pdf")

res_sg=results$`Single vs UT`
plot_top_dropouts(res_sg, filename_pdf="top40_singlet_tcr_vs_ut_dropouts.pdf")

# # individual genes
# -----------------------
# same size factors as in screen_calc()
counts_mat=counts_mat[rowSums(counts_mat)>30,] # slight different filtering
sf <- colSums(counts_mat) / max(colSums(counts_mat))

# divide each column (sample) by its size factor
norm_mat <- t(t(counts_mat) / sf)          # counts(eset, normalized=TRUE) equivalent

norm_dt <- as.data.table(norm_mat, keep.rownames = "guide")
mrow_n  <- melt(norm_dt, id.vars = "guide",
                variable.name = "variable", value.name = "norm")

# attach metadata (rows_dt already aligned to counts_mat via the match() earlier)
meta_cols <- c("guide","gene","gene_replicate","mutation_profile","guide_type",
               "gene_name","pep_id","gene_id","ref_key","ref_key_3p","ref_key_5p","is_fusion")
if ("barcode" %in% names(rows_dt)) setnames(rows_dt, "barcode", "guide")
mrow_n <- rows_dt[, ..meta_cols][mrow_n, on = "guide"]
fwrite(mrow_n, paste0(output_dir, "gene_profiles.txt"))

GENE_SHOW <- "TBX3"

rk <- mrow_n[gene_name == GENE_SHOW & !grepl("ref", mutation_profile), unique(ref_key)]
pd <- rbind(
  mrow_n[gene_name == GENE_SHOW & !grepl("ref", mutation_profile)],
  mrow_n[grepl("ref", mutation_profile) & ref_key %in% rk][
    unique(mrow_n[gene_name == GENE_SHOW & !grepl("ref", mutation_profile),
                  .(ref_key, panel = gene_id)]),
    on = "ref_key", allow.cartesian = TRUE],
  use.names = TRUE, fill = TRUE)
pd[is.na(panel), panel := gene_id]          # alt rows keep their own id

ggplot(pd, aes(variable, norm + 1)) +
  geom_point(aes(col = mutation_profile), size = 3, alpha = .5) +
  geom_line(aes(group = guide)) +
  facet_wrap(~ panel) + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_color_manual(values = palette_mutprof) +
  labs(x = NULL, y = "normalised counts + 1")
