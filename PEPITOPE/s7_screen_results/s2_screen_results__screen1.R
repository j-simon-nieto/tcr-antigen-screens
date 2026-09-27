# # REIVIST TO REMOVE UNNCESSARY VARIABLES FROM DESQE COMPARISONS> ALSO PLORBABLY REMOVE THE INCLUSION OF 2B in this way.
# # REPORT ON THE DECISIONS TAKEN THUS FAR


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

output_dir="PEPITOPE/s7_screen_results/output_screen_1/updated_screen_results/" 
# updates to the logfc thresholds in deseq2
# run on 18/Sept/26: fix fusion proteins
if(!dir.exists(output_dir)){
  dir.create(output_dir)}


# Grabbing TAA big library for barcodes

# 0/ Load files
# ================================================

# # results from QC
stats_dt=fread(paste0(output_dir, "../manual_qc_stats.txt"))
counts_dt=fread(paste0(output_dir, "../manual_qc_counts.txt"))
rows_dt=fread(paste0(output_dir, "../manual_qc_rows.txt"))

# # NEW! 18/09/26
# # # we have seen issues with the viability and killing capacity of UTs in this experiment
# # # the same controls were taken in 2b without said issues
# # # adding them as a "cleaner control" to see if it fixes the CMV issue (donor killing associated)
counts_2b_dt=fread(paste0(output_dir, "../../output_screen_2b/manual_qc_counts.txt"))
counts_mat_2b= as.matrix(as.data.frame(counts_2b_dt)[,3:ncol(counts_2b_dt)])
rownames(counts_mat_2b)=counts_2b_dt$guide
colnames(counts_mat_2b)=paste0("from2b__", colnames(counts_mat_2b))
counts_mat_2b=counts_mat_2b[, -ncol(counts_mat_2b)]

# 1/ Re-format and adapt DESEQ2
# ================================================

# create central counts and design
counts_mat= as.matrix(as.data.frame(counts_dt)[,3:ncol(counts_dt)])
rownames(counts_mat)=counts_dt$guide

# # join with screen 2b data
stopifnot(identical(rownames(counts_mat), rownames(counts_mat_2b)))
counts_mat=cbind(counts_mat,counts_mat_2b)

# # visualise comparisons
viz_counts=data.frame(counts_mat)
ggplot(viz_counts, aes(from2b__UT...C+1, UT...C+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
ggplot(viz_counts, aes(from2b__UT...B+1, UT...B+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
ggplot(viz_counts, aes(from2b__UT...A+1, UT...A+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
ggplot(viz_counts, aes(from2b__UT...A+1, Cluster.T...A+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
ggplot(viz_counts, aes(from2b__UT...A+1, Single.T...A+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
ggplot(viz_counts, aes(from2b__UT...A+1, X1D3...B+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")
# # seems like the 2b samples have lower counts among the most detected barcodes

ggplot(viz_counts, aes(B.only...A+1, from2b__B.only...A+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
  ggplot(viz_counts, aes(B.only...B+1, from2b__B.only...B+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
  ggplot(viz_counts, aes(B.only...C+1, from2b__B.only...C+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
  
  ggplot(viz_counts, aes(B.only...A+1, UT...A+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
  ggplot(viz_counts, aes(B.only...B+1, UT...B+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
  ggplot(viz_counts, aes(B.only...C+1, UT...C+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
  
  ggplot(viz_counts, aes(from2b__B.only...A+1, from2b__UT...A+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
  ggplot(viz_counts, aes(from2b__B.only...B+1, from2b__UT...B+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")+
  ggplot(viz_counts, aes(from2b__B.only...C+1, from2b__UT...C+1))+geom_point(alpha=.2,stroke=.2, size=.5)+scale_x_log10()+scale_y_log10()+geom_abline(slope=1, col="red")


design_dt=data.table(sample_id=colnames(counts_mat))
design_dt[, rep:= stringr::str_split_i(sample_id, " ", -1)]
design_dt[, origin:= stringr::str_split_i(sample_id, " ", 1)]

setnames(rows_dt, "guide", "barcode")
rows_dt <- rows_dt[match(rownames(counts_mat), rows_dt$barcode), ]
rows_dt[, gene:=gsub("^_","",gene)]
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
# -----------------------

results=screen_calc(dset, comparisons=list(c("UT", "B"),
                                           c("X1D3", "UT"),
                                           c("Single", "UT"),
                                           c("Cluster", "UT"),
                                           c("from2b__UT", "UT"),
                                           c("X1D3", "from2b__UT"),
                                           c("Single", "from2b__UT"),
                                           c("Cluster", "from2b__UT"),
                                           c("X1D3", "B"),
                                           c("Single", "B"),
                                           c("Cluster", "B"),
                                           c("Single", "X1D3"),
                                           c("Cluster", "X1D3")))

results_agg=rbindlist(lapply(results, setDT), idcol = "comparison")
fwrite(results_agg, paste0(output_dir, "screen1_results.txt"))

# Overview
# -----------------------

ggplot(results_agg[baseMean>100], aes(log2FoldChange, -log10(padj)))+geom_point(aes(col=mutation_profile))+
  ggrepel::geom_text_repel(data=results_agg[padj<1e-10],aes(label=gene_name))+
  facet_wrap(~comparison)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_vline(col="grey70", xintercept = c(-.1,.1), lty="dashed")
ggsave(paste0(output_dir, "volcano_plot_overview.pdf"), width = 12, height=8)

ggplot(results_agg[baseMean>0], aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_agg[padj<1e-5],aes(label=gene_name), size=2.2)+
  facet_wrap(~comparison)+scale_x_log10()+
  geom_point(shape=1, data=results_agg[padj<0.05], size=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_vline(col="grey70", xintercept = c(100), lty="dashed")
ggsave(paste0(output_dir, "MA_like_plot_overview.pdf"), width = 12, height=8)

ggplot(results_agg[baseMean>100 & (comparison%in% c("Cluster vs UT", "Single vs UT"))], aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_agg[baseMean>100 &padj<0.05& comparison%in% c("Cluster vs UT", "Single vs UT")],
                           aes(label=gene_id), size=2.5, max.overlaps = 15)+
  facet_wrap(~comparison)+scale_x_log10()+
  geom_point(shape=1, data=results_agg[baseMean>100 &padj<0.05& comparison%in% c("Cluster vs UT", "Single vs UT")], size=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))
ggsave(paste0(output_dir, "MA_like_plot_clusters_and_singlets_vs_UT.pdf"), width = 11, height=4)

ggplot(results_agg[grepl("^Clus|^Sin|^X1", comparison) & grepl("vs UT|from2b", comparison)], aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_agg[padj<1e-5 & grepl("^Clus|^Sin|^X1", comparison)  & grepl("vs UT|from2b", comparison)],
                           aes(label=gene_id), size=2.5, force = .3)+
  facet_wrap(~comparison, nrow=3)+scale_x_log10()+
  geom_point(shape=1, data=results_agg[padj<0.05 & grepl("^Clus|^Sin|^X1", comparison)  & grepl("vs UT|from2b", comparison)], size=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_hline(yintercept = c(0),lwd=1)
ggsave(paste0(output_dir, "MA_like_plot_clusters_and_singlets_and_1D3_vs_UT.pdf"), width = 8, height=8)

comp_db_sg_ut=merge(results$`Cluster vs UT`, results$`Single vs UT`, by=names(rows_dt)[!grepl(" ", names(rows_dt))], suffixes=c("_db", "_sg"))
ggplot(comp_db_sg_ut[baseMean_db>100 | baseMean_sg>100], aes(stat_sg, stat_db))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=comp_db_sg_ut[(baseMean_db>100 | baseMean_sg>100) & 
                                                (padj_db<0.05|padj_sg<0.05)],
                           aes(label=gene_id), size=3, max.overlaps = 20)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Cluster TCRs vs UT", x="Stat: Singlet TCRs vs UT")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")
ggsave(paste0(output_dir, "Clus_vs_UT__Singlet_vs_UT.pdf"), width = 5, height=5)

new_comp_db_sg_ut=merge(results$`Cluster vs from2b__UT`, results$`Single vs from2b__UT`, by=names(rows_dt)[!grepl(" ", names(rows_dt))], suffixes=c("_db", "_sg"))
ggplot(new_comp_db_sg_ut[baseMean_db>100 | baseMean_sg>100], aes(stat_sg, stat_db))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=new_comp_db_sg_ut[(baseMean_db>100 | baseMean_sg>100) & 
                                                (padj_db<0.05|padj_sg<0.05)],
                           aes(label=gene_id), size=2.4, max.overlaps = 30, force=.1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Cluster TCRs vs UT (s2B)", x="Stat: Singlet TCRs vs UT (s2B)")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")
ggsave(paste0(output_dir, "new_2b_controls_Clus_vs_UT__Singlet_vs_UT.pdf"), width = 5, height=5)

comp_db_ctrls=merge(results$`Cluster vs X1D3`, results$`Cluster vs UT`, by=names(rows_dt)[!grepl(" ", names(rows_dt))], suffixes=c("_1d3", "_ut"))
ggplot(comp_db_ctrls[baseMean_1d3>100 | baseMean_ut>100], aes(stat_ut, stat_1d3))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=comp_db_ctrls[(baseMean_1d3>100 | baseMean_ut>100) & 
                                                (padj_1d3<0.05|padj_ut<0.05)],
                           aes(label=gene_id), size=3, max.overlaps = 20)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Cluster TCRs vs 1D3", x="Stat: Cluster TCRs vs UT")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")
ggsave(paste0(output_dir, "Clus_vs_UT__Cluster_vs_1D3.pdf"), width = 5, height=5)

comp_sg_ctrls=merge(results$`Single vs X1D3`, results$`Single vs UT`, 
                    by=names(rows_dt)[!grepl(" ", names(rows_dt))], suffixes=c("_1d3", "_ut"))
ggplot(comp_sg_ctrls[baseMean_1d3>100 | baseMean_ut>100], aes(stat_ut, stat_1d3))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=comp_sg_ctrls[(baseMean_1d3>100 | baseMean_ut>100) & 
                                                (padj_1d3<0.05|padj_ut<0.05)],
                           aes(label=gene_id), size=3, max.overlaps = 20)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Single TCRs vs 1D3", x="Stat: Single TCRs vs UT")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")
ggsave(paste0(output_dir, "Single_vs_UT__Singlet_vs_1D3.pdf"), width = 5, height=5)

comp_sg_ctrls_2b=merge(results$`Single vs from2b__UT`, results$`Single vs UT`, 
                    by=names(rows_dt)[!grepl(" ", names(rows_dt))], suffixes=c("_from_2b_ut", "_ut"))
ggplot(comp_sg_ctrls_2b[baseMean_from_2b_ut>100 | baseMean_ut>100], aes(stat_ut,stat_from_2b_ut))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=comp_sg_ctrls_2b[(baseMean_from_2b_ut>100 | baseMean_ut>100) & 
                                                (padj_from_2b_ut<0.05|padj_ut<0.05)],
                           aes(label=gene_id), size=2.4, max.overlaps = 20)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Single TCRs vs UT (s2B)", x="Stat: Single TCRs vs UT")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")
ggsave(paste0(output_dir, "Single_vs_UT__Singlet_vs_UT_from_2b.pdf"), width = 5, height=5)

comp_db_ctrls_2b=merge(results$`Cluster vs from2b__UT`, results$`Cluster vs UT`, 
                       by=names(rows_dt)[!grepl(" ", names(rows_dt))], suffixes=c("_from_2b_ut", "_ut"))
ggplot(comp_db_ctrls_2b[baseMean_from_2b_ut>100 | baseMean_ut>100], aes(stat_ut,stat_from_2b_ut))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=comp_db_ctrls_2b[(baseMean_from_2b_ut>100 | baseMean_ut>100) & 
                                                (padj_from_2b_ut<0.05|padj_ut<0.05)],
                           aes(label=gene_id), size=2.4, max.overlaps = 20)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Cluster TCRs vs UT (s2B)", x="Stat: Cluster TCRs vs UT")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")
ggsave(paste0(output_dir, "Cluster_vs_UT__Clustert_vs_UT_from_2b.pdf"), width = 5, height=5)

# # Analysis of top hits:
# -----------------------
plot_top_dropouts=
  function(res_dt, top_n=40, filename_pdf){

  res_dt=res_dt[baseMean>30,]
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

res_clus=results$`Cluster vs X1D3`
plot_top_dropouts(res_clus, filename_pdf="top40_cluster_tcr_vs_1d3_dropouts.pdf")

res_sg=results$`Single vs X1D3`
plot_top_dropouts(res_sg, filename_pdf="top40_singlet_tcr_vs_1d3_dropouts.pdf")

# # individual genes
# -----------------------
# same size factors as in screen_calc()
counts_mat=counts_mat[rowSums(counts_mat)>30,]# slight different filtering
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



# DIAGNOSIS: we see population shifts, we want to explore where are they coming from

design_dt[, grp := factor(paste(fifelse(grepl("^from2b__", sample_id), "s2b", "s1"),
                                sub("^from2b__", "", origin), sep="."))]

dds <- DESeqDataSetFromMatrix(counts_mat, design_dt, ~ 0 + grp, rowData=rows_dt)
dds <- dds[matrixStats::rowMaxs(counts(dds)) >= 30, ]
sizeFactors(dds) <- colSums(counts(dds)) / max(colSums(counts(dds)))
dds <- DESeq(dds, fitType="local")

did <- function(a, b, c, d) {                      # (a-b) - (c-d)
  v <- setNames(numeric(length(resultsNames(dds))), resultsNames(dds))
  v[paste0("grp", c(a, b, c, d))] <- c(1, -1, -1, 1); v
}

chk <- function(r) as.data.table(r)[!is.na(padj), .(n=.N, nsig=sum(padj < 0.05),
                                                    medLFC=round(median(log2FoldChange), 3),
                                                    trend=round(cor(log10(baseMean), log2FoldChange), 3))]

rbind(placebo = chk(results(dds, contrast=did("s1.1D3","s1.B","s2b.1D3","s2b.B"))),
      ut_delta = chk(results(dds, contrast=did("s1.UT","s1.B","s2b.UT","s2b.B"))),
      idcol="test")

chk(results(dds, contrast = did("s1.1D3","s2b.1D3","s1.UT","s2b.UT")))
