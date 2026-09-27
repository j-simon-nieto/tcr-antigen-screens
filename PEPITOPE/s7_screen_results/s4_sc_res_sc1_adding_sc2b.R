# SEE RESULTS IN TCR_antigen_screen_hybdrization
library(ggrepel)
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

output_dir="PEPITOPE/s7_screen_results/output_screen_1/hybridized_screen_results/" 
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

# # remove all variables that we are not using for comparisons
design_dt[, clean_origin:= gsub("from2b__","",origin)]
design_dt_f=design_dt[!clean_origin %in% c("unmatched")] # remove Bs eventually -- keeping them as controls now
design_dt_f[, screen:= ifelse(grepl("from2b", origin), "s2b", "s1")]
counts_mat_f=counts_mat[, colnames(counts_mat) %in% design_dt_f$sample_id]

dset <- DESeqDataSetFromMatrix(counts_mat_f, design_dt_f, ~ clean_origin, rowData=rows_dt)

# # Deseq function
# Deseq function
screen_calc = function(dset, comparisons, min_count=30, two_screens=T) {
  if (length(unique(dset$patient)) > 1)
    stop("The 'patient' column can not span more than one value")
  if (is.numeric(dset$rep))
    dset$rep = factor(dset$rep)
  
  # eset = DESeq2::DESeqDataSet(dset, ~ rep + origin) # removing because the replicates are not paired!!!
  
  if(two_screens==F){
    eset = DESeq2::DESeqDataSet(dset, ~ clean_origin)
  }else{
  eset = DESeq2::DESeqDataSet(dset, ~ screen + clean_origin)}
  
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
    DESeq2::results(mod, contrast=c("clean_origin", comp)) |>
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
results=screen_calc(dset, comparisons=list(c("1D3", "UT"),
                                           c("Single", "UT"),
                                           c("Cluster", "UT"),
                                           c("Single", "1D3"),
                                           c("Cluster", "1D3"),
                                           c("UT", "B")))

results_agg=rbindlist(lapply(results, setDT), idcol = "comparison")
fwrite(results_agg, paste0(output_dir, "hyb_screen1_results.txt"))

# Viz DESEQ results
# -----------------------
ggplot(results_agg[baseMean>100], aes(log2FoldChange, -log10(padj)))+geom_point(aes(col=mutation_profile))+
  ggrepel::geom_text_repel(data=results_agg[padj<1e-10],aes(label=gene_name))+
  facet_wrap(~comparison)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))#+
  #geom_vline(col="grey70", xintercept = c(-.1,.1), lty="dashed")
#ggsave(paste0(output_dir, "volcano_plot_overview.pdf"), width = 12, height=8)

ggplot(results_agg[comparison=="Cluster vs UT"], aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_agg[comparison=="Cluster vs UT" &padj<1e-5],aes(label=gene_name), size=2.2)+
  facet_wrap(~comparison)+scale_x_log10()+
  geom_point(shape=1, data=results_agg[comparison=="Cluster vs UT" & padj<0.05], size=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_vline(col="grey70", xintercept = c(100), lty="dashed")
#ggsave(paste0(output_dir, "MA_like_plot_overview.pdf"), width = 12, height=8)

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

comp_db_ctrls=merge(results$`Cluster vs 1D3`, results$`Cluster vs UT`, by=names(rows_dt)[!grepl(" - ", names(rows_dt))], suffixes=c("_1d3", "_ut"))
ggplot(comp_db_ctrls[baseMean_1d3>100 | baseMean_ut>100], aes(stat_ut, stat_1d3))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=comp_db_ctrls[(baseMean_1d3>100 | baseMean_ut>100) & 
                                                (padj_1d3<0.05|padj_ut<0.05)],
                           aes(label=gene_id), size=4, max.overlaps = 20)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Cluster TCRs vs 1D3", x="Stat: Cluster TCRs vs UT")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")


# PREVIOUS ITERATIONS:
# ---------------------

# # ONLY USING S1
dess1=design_dt_f[screen=="s1"]
cms1=counts_mat_f[, colnames(counts_mat_f) %in% dess1$sample_id]
dset_s1 <- DESeqDataSetFromMatrix(cms1, dess1, ~ clean_origin, rowData=rows_dt)

results_s1=screen_calc(dset_s1, comparisons=list(c("1D3", "UT"),
                                           c("Single", "UT"),
                                           c("Cluster", "UT"),
                                           c("Single", "1D3"),
                                           c("Cluster", "1D3"),
                                           c("UT", "B")), two_screens = F)

results_s1_agg=rbindlist(lapply(results_s1, setDT), idcol = "comparison")

s1_comp_db_ctrls=merge(results_s1$`Cluster vs 1D3`, results_s1$`Cluster vs UT`, by=names(rows_dt)[!grepl(" - ", names(rows_dt))], suffixes=c("_1d3", "_ut"))
ggplot(s1_comp_db_ctrls[baseMean_1d3>100 | baseMean_ut>100], aes(stat_ut, stat_1d3))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=s1_comp_db_ctrls[(baseMean_1d3>100 | baseMean_ut>100) & 
                                                (padj_1d3<0.05|padj_ut<0.05)],
                           aes(label=gene_id), size=4, max.overlaps = 30)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Cluster TCRs vs 1D3", x="Stat: Cluster TCRs vs UT")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")

ggplot(results_s1_agg[baseMean>100], aes(log2FoldChange, -log10(padj)))+geom_point(aes(col=mutation_profile))+
  ggrepel::geom_text_repel(data=results_s1_agg[padj<1e-10],aes(label=gene_name))+
  facet_wrap(~comparison)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))

# # TRASNPOSING THE S2B CONTROL
dess2b=copy(design_dt_f)
dess2b[, clean_origin:= origin]
cms2b=copy(counts_mat_f)
dset_s2b <- DESeqDataSetFromMatrix(cms2b, dess2b, ~ clean_origin, rowData=rows_dt)

results_s2b=screen_calc(dset_s2b, comparisons=list(c("1D3", "UT"),
                                                 c("Single", "UT"),
                                                 c("Cluster", "UT"),
                                                 c("1D3", "from2b__UT"),
                                                 c("Single", "from2b__UT"),
                                                 c("Cluster", "from2b__UT"),
                                                 c("Single", "1D3"),
                                                 c("Cluster", "1D3"),
                                                 c("UT", "B")), two_screens = F)

results_s2b_agg=rbindlist(lapply(results_s2b, setDT), idcol = "comparison")

s2b_comp_db_ctrls=merge(results_s2b$`Cluster vs from2b__UT`, results_s2b$`Cluster vs UT`, by=names(rows_dt)[!grepl(" - ", names(rows_dt))], suffixes=c("_s2b_ut", "_s1_ut"))
ggplot(s2b_comp_db_ctrls[baseMean_s1_ut>100 | baseMean_s2b_ut>100], aes(stat_s1_ut, stat_s2b_ut))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=s2b_comp_db_ctrls[(baseMean_s1_ut>100 | baseMean_s2b_ut>100) & 
                                                   grepl("FLU|CMV",gene_name)],
                           aes(label=gene_id), size=3, max.overlaps = 30)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Cluster TCRs vs UT_s2b", x="Stat: Cluster TCRs vs UT_s1")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")

ggplot(results_s2b_agg[comparison %in% c("Cluster vs from2b__UT", "Cluster vs UT")],
       aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_s2b_agg[comparison %in% c("Cluster vs from2b__UT", "Cluster vs UT")][padj<1e-5],
                           aes(label=gene_name), size=2.2)+
  facet_wrap(~comparison)+scale_x_log10()+
  geom_point(shape=1, data=results_s2b_agg[comparison %in% c("Cluster vs from2b__UT", "Cluster vs UT")][padj<0.05], size=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_vline(col="grey70", xintercept = c(100), lty="dashed")

# counts distribution
cms_2b_dt=data.table(cms2b)

cms_2b_dt[, mean_UT_s1:= rowMeans(.SD), .SDcols=c("UT - A", "UT - B", "UT - C")]
cms_2b_dt[, mean_UT_s2b:= rowMeans(.SD), .SDcols=c("from2b__UT - A", "from2b__UT - B", "from2b__UT - C")]

ggplot(cms_2b_dt, aes(mean_UT_s1,mean_UT_s2b))+
  geom_point(alpha=.2,stroke=.3, size=1)+
  scale_x_log10()+scale_y_log10()+
  geom_abline(slope=1, col="red")+
  ggtitle("Barcode counts across screens\n(UT condition)")+

ggplot(cms_2b_dt, aes(mean_UT_s1,mean_UT_s2b))+
  geom_point(alpha=.2,stroke=.3, size=1)+
  #scale_x_log10()+scale_y_log10()+
  geom_abline(slope=1, col="red")+
  ggtitle("Barcode counts across screens\n(UT condition)")

# comparing the solutions
comp_old_new=merge(results_s1$`Cluster vs UT`, results$`Cluster vs UT`, by=names(rows_dt)[!grepl(" - ", names(rows_dt))], suffixes=c("_old", "_new"))
problematic=c("CMV_pp65-27", "CMV_pp65-26", "FLU_NP-27")
ggplot(comp_old_new[baseMean_old>100 & baseMean_new>100],
       aes(log2FoldChange_old, log2FoldChange_new))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_text(data=comp_old_new[gene_id %in% problematic], aes(label=gene_id))

ggplot(comp_old_new[baseMean_old>100 & baseMean_new>100],
       aes(stat_old, stat_new))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_text_repel(data=comp_old_new[gene_id %in% problematic], aes(label=gene_id))


# A NEW DIRECTION: can we remove the old UT and only use Bcells and 1D3 as a batch effect alignment fraction
design_dt_f_2=design_dt_f[!(screen=="s1" & origin=="UT")]
counts_mat_f_2=counts_mat_f[, colnames(counts_mat_f) %in% design_dt_f_2$sample_id]

dset_hyb <- DESeqDataSetFromMatrix(counts_mat_f_2, design_dt_f_2, ~ clean_origin, rowData=rows_dt)

results_hyb=screen_calc(dset_hyb, comparisons=list(c("1D3", "UT"),
                                                 c("Single", "UT"),
                                                 c("Cluster", "UT"),
                                                 c("Single", "1D3"),
                                                 c("Cluster", "1D3"),
                                                 c("UT", "B")))

results_hyb_agg=rbindlist(lapply(results_hyb, setDT), idcol = "comparison")

# Viz DESEQ results
# -----------------------
setDT(results_hyb_agg)
ggplot(results_hyb_agg[baseMean>100], aes(log2FoldChange, -log10(padj)))+geom_point(aes(col=mutation_profile))+
  ggrepel::geom_text_repel(data=results_hyb_agg[padj<1e-10],aes(label=gene_name))+
  facet_wrap(~comparison)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))#+
#geom_vline(col="grey70", xintercept = c(-.1,.1), lty="dashed")
#ggsave(paste0(output_dir, "volcano_plot_overview.pdf"), width = 12, height=8)

ggplot(results_hyb_agg[comparison=="Cluster vs UT"], aes(baseMean, log2FoldChange))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  ggrepel::geom_text_repel(data=results_hyb_agg[comparison=="Cluster vs UT" &padj<1e-5],aes(label=gene_name), size=2.2)+
  facet_wrap(~comparison)+scale_x_log10()+
  geom_point(shape=1, data=results_hyb_agg[comparison=="Cluster vs UT" & padj<0.05], size=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_vline(col="grey70", xintercept = c(100), lty="dashed")
#ggsave(paste0(output_dir, "MA_like_plot_overview.pdf"), width = 12, height=8)

comp_db_sg_ut_hyb=merge(results_hyb$`Cluster vs UT`, results_hyb$`Single vs UT`,
                    by=names(rows_dt)[!grepl(" ", names(rows_dt))], suffixes=c("_db", "_sg"))
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

comp_db_ctrls=merge(results_hyb$`Cluster vs 1D3`, results_hyb$`Cluster vs UT`,
                    by=names(rows_dt)[!grepl(" - ", names(rows_dt))], suffixes=c("_1d3", "_ut"))
ggplot(comp_db_ctrls[baseMean_1d3>100 | baseMean_ut>100], aes(stat_ut, stat_1d3))+
  geom_point(aes(col=mutation_profile), size = 2, alpha=1)+
  ggrepel::geom_text_repel(data=comp_db_ctrls[(baseMean_1d3>100 | baseMean_ut>100) & 
                                                (padj_1d3<0.05|padj_ut<0.05)],
                           aes(label=gene_id), size=4, max.overlaps = 20)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  labs(y="Stat: Cluster TCRs vs 1D3", x="Stat: Cluster TCRs vs UT")+
  theme(legend.position = "none", axis.title = element_text(size=14))+
  geom_hline(yintercept = 0, lty="dashed", col="grey")+
  geom_vline(xintercept = 0, lty="dashed", col="grey")


comp_old_new=merge(results_s1$`Cluster vs UT`, results_hyb$`Cluster vs UT`, by=names(rows_dt)[!grepl(" - ", names(rows_dt))], suffixes=c("_old", "_new"))
problematic=c("CMV_pp65-27", "CMV_pp65-26", "FLU_NP-27")
ggplot(comp_old_new[baseMean_old>100 & baseMean_new>100],
       aes(log2FoldChange_old, log2FoldChange_new))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_text_repel(data=comp_old_new[gene_id %in% problematic], aes(label=gene_id))

ggplot(comp_old_new[baseMean_old>100 & baseMean_new>100],
       aes(stat_old, stat_new))+
  geom_point(aes(col=mutation_profile), size = 1, alpha=1)+
  scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold"))+
  geom_text_repel(data=comp_old_new[gene_id %in% problematic], aes(label=gene_id))+
  geom_abline(slope=1)


# VISUALISE HOW THE BATCH INTEGRATION SHITS DEPENDING ON THE OVERLAP
# ============================================================
# Effect of keeping vs dropping B as a batch-alignment anchor
# ============================================================

# --- model 1: B kept (= current results_hyb) --------------------------------
res_wB <- setDT(copy(results_hyb$`Cluster vs UT`))

# --- model 2: B dropped (1D3 as sole anchor) --------------------------------
design_dt_f_3 <- design_dt_f_2[clean_origin != "B"]
counts_mat_f_3 <- counts_mat_f[, colnames(counts_mat_f) %in% design_dt_f_3$sample_id]
dset_noB <- DESeqDataSetFromMatrix(counts_mat_f_3, design_dt_f_3, ~ clean_origin, rowData = rows_dt)
res_nB <- setDT(screen_calc(dset_noB, comparisons = c("Cluster", "UT")))

# --- model 3: unmerged origins, to measure each anchor's screen offset ------
design_dt_f_4 <- copy(design_dt_f_2)[, clean_origin := origin]
counts_mat_f_4 <- counts_mat_f[, colnames(counts_mat_f) %in% design_dt_f_4$sample_id]
dset_anch <- DESeqDataSetFromMatrix(counts_mat_f_4, design_dt_f_4, ~ clean_origin, rowData = rows_dt)
res_anch <- screen_calc(dset_anch,
                        comparisons = list(c("from2b__B", "B"), c("from2b__1D3", "1D3")),
                        two_screens = FALSE)

# --- assemble ---------------------------------------------------------------
cmp <- merge(res_wB[, .(barcode, gene_id, mutation_profile, baseMean,
                        lfc_wB = log2FoldChange, se_wB = lfcSE, padj_wB = padj)],
             res_nB[, .(barcode, lfc_nB = log2FoldChange, se_nB = lfcSE, padj_nB = padj)],
             by = "barcode")
cmp <- merge(cmp, setDT(res_anch$`from2b__B vs B`)[, .(barcode, dB = log2FoldChange)], by = "barcode")
cmp <- merge(cmp, setDT(res_anch$`from2b__1D3 vs 1D3`)[, .(barcode, d1 = log2FoldChange)], by = "barcode")
cmp[, `:=`(delta = lfc_wB - lfc_nB, donor = gene_id %in% problematic)]

message(sprintf("barcodes: wB=%d  nB=%d  merged=%d", nrow(res_wB), nrow(res_nB), nrow(cmp)))
message(sprintf("delta~baseMean trend = %.3f | median SE ratio (wB/nB) = %.3f | nsig %d vs %d",
                cmp[baseMean > 30, cor(log10(baseMean), delta, use = "complete")],
                cmp[, median(se_wB / se_nB, na.rm = TRUE)],
                cmp[, sum(padj_wB < 0.05, na.rm = TRUE)], cmp[, sum(padj_nB < 0.05, na.rm = TRUE)]))

pal <- c(alt = "red3", ref = "grey", "TAA" = "gold")

p1 <- ggplot(cmp[baseMean > 100], aes(lfc_nB, lfc_wB)) +
  geom_point(aes(col = mutation_profile), size = 1, alpha = .5) +
  geom_abline(col = "grey40") + coord_equal() +
  geom_text_repel(data = cmp[baseMean > 100 & donor], aes(label = gene_id), size = 3) +
  scale_color_manual(values = pal) + theme(legend.position = "none") +
  labs(x = "Cluster vs UT  (1D3 anchor only)", y = "Cluster vs UT  (B + 1D3 anchors)",
       title = "Do the two models agree?")

p2 <- ggplot(cmp[baseMean > 100], aes(log10(baseMean), delta)) +
  geom_point(alpha = .12, size = .6) +
  geom_hline(yintercept = 0, col = "grey50") +
  geom_smooth(method = "loess", col = "red", se = FALSE) +
  labs(x = "log10 baseMean", y = "delta log2FC  (B kept - B dropped)",
       title = "Does B inject the abundance tilt?")

p3 <- ggplot(cmp, aes(0.5 * (dB - d1), delta)) +
  geom_point(alpha = .15, size = .6) +
  geom_abline(col = "red") + coord_equal() +
  labs(x = "1/2 * (dB - d1D3)", y = "observed delta",
       title = "Mechanism: half the anchor disagreement")

p4 <- ggplot(melt(cmp[donor == TRUE & baseMean > 100,
                      .(gene_id, `1D3 only` = lfc_nB, `B + 1D3` = lfc_wB)], id.vars = "gene_id"),
             aes(value, reorder(gene_id, value))) +
  geom_line(aes(group = gene_id), col = "grey70") +
  geom_point(aes(col = variable), size = 3) +
  geom_vline(xintercept = 0, lty = "dashed", col = "grey") +
  labs(x = "log2 fold change", y = NULL, col = NULL, title = "Donor-reactive barcodes")

cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv", axis = "lrtb")
#ggsave(paste0(output_dir, "anchor_B_vs_1D3_comparison.pdf"), width = 11, height = 9)
