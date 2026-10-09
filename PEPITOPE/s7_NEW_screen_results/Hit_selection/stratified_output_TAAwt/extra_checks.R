library(data.table)

# Grab barcodes
#----------------------------------
barcode_ev_p5=fread("PEPITOPE/s7_NEW_screen_results/Hit_selection/stratified_output_TAAwt/P5/04_barcode_evidence.txt")

barcode_ev_p5[gene_id=="PROSER1_P231S"][comparison!="UT vs B"]
barcode_ev_p5[gene_id=="ERVE-1.2-2"][comparison!="UT vs B"]

screen_p5=fread("PEPITOPE/s7_NEW_screen_results/output_screen_2a/updated_screen_results/screen2a_results.txt")
screen_p5[gene_id=="ERVE-1.2-1"][comparison %in% c("Cluster vs UT","Single vs UT")]
screen_p5[gene_id=="P53-8"][comparison %in% c("Cluster vs UT","Single vs UT")]

# easy visualisation of cluster/singlet
#----------------------------------
screen_p1=fread("PEPITOPE/s7_NEW_screen_results/output_screen_1/updated_screen_results/screen1_results.txt")
screen_p5=fread("PEPITOPE/s7_NEW_screen_results/output_screen_2a/updated_screen_results/screen2a_results.txt")
selection=fread('PEPITOPE/s7_NEW_screen_results/Hit_selection/stratified_output_TAAwt/Combined_selection.csv', skip = 0, header = T)
source("PEPITOPE/s7_NEW_screen_results/palette_screens.R")

# (1) define the helper your cell_fun calls
.sig_stars <- function(p) ifelse(is.na(p), "",
                                 ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

make_heat_compare_comps=function(screen, pat){
  tot_res=screen[,.(pat=pat,barcode, mutation_profile, comparison, log2FoldChange, padj, gene_id)]
  selection=selection[Patient==pat,]
  sel_barcodes=tot_res[barcode %in% selection[`Clone?`=="yes"]$`barcode JSN`]
  uniqueN(sel_barcodes$gene_id) #24
  
  f_res=tot_res[gene_id %in% sel_barcodes$gene_id][comparison %in% c("Cluster vs UT","Single vs UT")][pat==pat]
  f_res[, selected_barcode:= ifelse(barcode %in% sel_barcodes$barcode, "Selected", "Other")]
  f_res[, selected_barcode:= make.unique(selected_barcode), by=c("gene_id", "comparison")]
  
  library(ggplot2)
  heatmap_l2fc=dcast(f_res, gene_id+mutation_profile~comparison+selected_barcode, value.var = "log2FoldChange")
  heatmap_padj=dcast(f_res, gene_id~comparison+selected_barcode, value.var = "padj")
  
  library(ComplexHeatmap)
  library(circlize)
  pal=colorRamp2(breaks=c(-3, -1, 0, 1, 3), colors = c("blue4", "royalblue1","white","red1", "red4"))
  
  row_anno=rowAnnotation(df=data.frame(mutation_profile=heatmap_l2fc$mutation_profile),
                         col=list(mutation_profile=palette_mutprof))
  gene_names=heatmap_l2fc$gene_id
  heatmap_l2fc=heatmap_l2fc[, 3:ncol(heatmap_l2fc)]
  
  top_anno= HeatmapAnnotation(df=data.frame(comparison=stringr::str_split_i(names(heatmap_l2fc),"_",1),
                                            barcode=stringr::str_split_i(names(heatmap_l2fc),"_",2)),
                              col = list(comparison=c("Cluster vs UT"="red4", "Single vs UT"="steelblue1"),
                                         barcode=c("Selected"="black", "Other"="grey")))
  
  heatmap_l2fc=data.frame(heatmap_l2fc, row.names = gene_names)
  
  # (2) the two matrices cell_fun indexes. padj value-columns come out of dcast in
  #     the SAME order as the l2fc value-columns, and rows are both sorted by
  #     gene_id, so they align positionally (cell_fun uses original-matrix i,j).
  mat  <- as.matrix(heatmap_l2fc)
  pmat <- as.matrix(heatmap_padj[, -1])
  
  cell_fun <- function(j, i, x, y, w, h, fill) {
    lab <- .sig_stars(pmat[i, j])
    if (nzchar(lab)) {
      v <- mat[i, j]
      col <- if (!is.na(v) && abs(v) > 0.35) "white" else "black"
      grid.text(lab, x, y, gp = gpar(fontsize = 9, col = col))
    }
  }
  
  sig_lgd <- Legend(title = "padj", type = "points", pch = NA,
                    labels = c("* <0.05", "** <0.01", "*** <0.001", "**** <0.00001"),
                    legend_gp = gpar(col = NA))
  
  # (3) pass cell_fun, and (4) draw with the significance legend
  ht <- Heatmap(mat, top_annotation = top_anno, col=pal, name="L2FC",
                left_annotation = row_anno, cell_fun = cell_fun)
  draw(ht, annotation_legend_list = list(sig_lgd), merge_legend = TRUE)
}

output_dir="PEPITOPE/s7_NEW_screen_results/Hit_selection/stratified_output_TAAwt/extra_checks/"
pdf(width = 5, height = 5, paste0(output_dir, "cluster_or_singlet_patient1.pdf"))
make_heat_compare_comps(screen_p1, "P1")
dev.off()

pdf(width = 6.5, height = 5, paste0(output_dir, "cluster_or_singlet_patient5.pdf"))
make_heat_compare_comps(screen_p5, "P5")
dev.off()

