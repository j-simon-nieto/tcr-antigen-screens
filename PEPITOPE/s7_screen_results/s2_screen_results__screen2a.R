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

output_dir="PEPITOPE/s7_screen_results/output_screen_2a/screen_results/"
if(!dir.exists(output_dir)){
  dir.create(output_dir)}

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

setnames(rows_dt, "guide", "barcode")
rows_dt <- rows_dt[match(rownames(counts_mat), rows_dt$barcode), ]
rows_dt[, gene_replicate:=stringr::str_split_i(gene, "_", -1)]
rows_dt[, mutation_profile:=stringr::str_split_i(gene, "_", -2)]
rows_dt[, gene_id:=stringr::str_split_i(gene, "_", 1)]
rows_dt[, gene_name:=gsub("^_", "", stringr::str_split_i(gene_id, "\\-", 1))]

rows_dt[, gene_id:= ifelse(mutation_profile=="TAA",
                           stringr::str_split_i(gene, "_TAA_TAA", 1), 
                           gene_id)]
rows_dt[, gene_id:= ifelse(gene_id=="",
                           gsub("^_","",stringr::str_split_i(gene, "--",1)),
                           gene_id)]

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
  DESeq2::sizeFactors(eset) = colSums(assay(dset)) / max(colSums(assay(dset)))
  mod = DESeq2::DESeq(eset, fitType="local")
  
  get_result = function(comp) {
    DESeq2::results(mod, contrast=c("origin", comp)) |>
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

ggplot(results_agg[baseMean>100], aes(log2FoldChange, -log10(padj)))+geom_point(aes(col=mutation_profile))+
  ggrepel::geom_text_repel(data=results_agg[padj<1e-10],aes(label=gene_name))+
  facet_wrap(~comparison)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))
ggsave(paste0(output_dir, "volcano_plot_overview.pdf"), width = 12, height=8)

ggplot(results_agg[baseMean>100], aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_agg[padj<1e-5],aes(label=gene_name), size=2.2)+
  facet_wrap(~comparison)+scale_x_log10()+
  geom_point(shape=1, data=results_agg[padj<0.05], size=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))
ggsave(paste0(output_dir, "MA_like_plot_overview.pdf"), width = 12, height=8)

ggplot(results_agg[baseMean>100 & grepl("vs UT", comparison)], aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_agg[padj<1e-5 & grepl("vs UT", comparison)],aes(label=gene_id), size=2.2)+
  facet_wrap(~comparison)+scale_x_log10()+
  geom_point(shape=1, data=results_agg[padj<0.05 & grepl("vs UT", comparison)], size=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_hline(yintercept = 0)

comp_db_sg_ut=merge(results$`Cluster vs UT`, results$`Single vs UT`, by=names(rows_dt)[!grepl(" ", names(rows_dt))], suffixes=c("_db", "_sg"))

ggplot(comp_db_sg_ut[baseMean_db>100 | baseMean_sg>100], aes(stat_sg, stat_db))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=comp_db_sg_ut[(baseMean_db>100 | baseMean_sg>100) & 
                                                (padj_db<0.05|padj_sg<0.05)],
                           aes(label=gene_name), size=3, max.overlaps = 10)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Cluster TCRs vs 1D3", x="Stat: Singlet TCRs vs 1D3")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")
ggsave(paste0(output_dir, "Clus_vs_UT__Singlet_vs_UT.pdf"), width = 5, height=5)

# # Analysis of top hits:
plot_top_dropouts=
  function(res_dt, top_n=40, filename_pdf){
    res_dt[, gene_id:= ifelse(mutation_profile=="TAA",
                              stringr::str_split_i(gene, "_TAA_TAA", 1), 
                              gene_id)]
    
    
    res_dt[, gene_id:= ifelse(gene_id=="",
                              gsub("^_","",stringr::str_split_i(gene, "--",1)),
                              gene_id)]
    res_dt=res_dt[baseMean>100,]
    setorder(res_dt, stat)
    top30=unique(res_dt$gene_id)[1:top_n]
    
    res_dt[, min_stat:= min(stat), by="gene_id"]
    sd_x=sd(res_dt$stat)
    mean_x=mean(res_dt$stat)
    
    ggplot(res_dt[gene_id %in% top30], 
           aes(reorder(gene_id, -min_stat), stat, col=mutation_profile))+
      geom_point(size=0.8)+coord_flip()+
      ggbeeswarm::geom_beeswarm(cex=.15)+
      geom_hline(yintercept = c(mean_x-sd_x, mean_x, mean_x+sd_x), 
                 lty=c("dashed", "solid","dashed"),
                 lwd=c(.3,.5,.3))+
      scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
      geom_point(shape=1, data=res_dt[gene_id %in% top30 & padj<0.05], size=3, col="black")+
      theme(legend.position = "none", axis.title = element_text(size=14))+
      scale_y_reverse()+xlab("")+ylab("Drop-out confidence (T-stat)")
    ggsave(paste0(output_dir, filename_pdf), width = 6, height=6)
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
meta_cols <- c("guide","gene","gene_replicate","mutation_profile",
               "gene_id","gene_name","guide_type")
setnames(rows_dt, "barcode", "guide")
mrow_n <- rows_dt[, ..meta_cols][mrow_n, on = "guide"]

fwrite(mrow_n, paste0(output_dir,"gene_profiles.txt"))

ggplot(mrow_n[gene_id == "DIDO1"], aes(variable, log10(norm + 1))) +
  geom_point(aes(col = mutation_profile), size = 4, alpha = .5) +
  geom_line(aes(group = guide)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  scale_color_manual(values = c(alt = "red3", ref = "grey", "TAA" = "gold"))
