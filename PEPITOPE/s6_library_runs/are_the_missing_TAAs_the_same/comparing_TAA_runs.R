library(data.table)
library(ggplot2)

# 0/ load data
# ==================
march_big=fread("PEPITOPE/s6_library_runs/first_run_march_26/output/report_source_data_pat_lib_plus_bigTAA.txt")
march_small=fread("PEPITOPE/s6_library_runs/first_run_march_26/output/report_source_data_pat_lib_plus_smallTAA.txt")

may_1_big=fread("PEPITOPE/s6_library_runs/second_run_may_26/output_sample_1/report_source_data_pat_lib_plus_bigTAA.txt")
may_1_small=fread("PEPITOPE/s6_library_runs/second_run_may_26/output_sample_1/report_source_data_pat_lib_plus_smallTAA.txt")

may_2_big=fread("PEPITOPE/s6_library_runs/second_run_may_26/output_sample_2/report_source_data_pat_lib_plus_bigTAA.txt")
may_2_small=fread("PEPITOPE/s6_library_runs/second_run_may_26/output_sample_2/report_source_data_pat_lib_plus_smallTAA.txt")

# # get big
big_taa=rbindlist(list("old_prep_original"=march_big,
               "old_prep_concentrated"=may_1_big,
               "new_prep"=may_2_big), fill=T, idcol = "preparation")
big_taa=big_taa[patient_id=="TAA_big" & detected_in=="TAA_big"]

# # get small
small_taa=rbindlist(list("old_prep_original"=march_small,
                       "old_prep_concentrated"=may_1_small,
                       "new_prep"=may_2_small), fill=T, idcol = "preparation")
small_taa=small_taa[patient_id=="TAA_small" & detected_in=="TAA_small"]

# # dictionary
dict_taa_big=fread("PEPITOPE/s6_library_runs/shared_by_JV/large_TAA.csv")

big_taa=merge(big_taa[, gene_name:=NULL], 
              dict_taa_big[, .(barcode=BARCODE,gene_name=`GENE NAME 2`)])

big_taa[, gene_name:= sub("-[^-]+$", "", gene_name)]


# # combine
libr_dt=rbind(small_taa[, library:="small_TAA"],big_taa[, library:="big_TAA"], fill=T) 
libr_dt[, count_lvl:= cut(counts, breaks = c(0,1,10,100,300,500,1000,Inf), include.lowest = T, ordered_result = T)]
fix_libr_dt=unique(libr_dt[, rank:=NULL][, rank_norm:=NULL])

# 1/ Plot data
# ==================

ggplot(fix_libr_dt, aes(preparation, counts+1))+geom_boxplot()+scale_y_log10()+
  theme_classic()+facet_wrap(~library)+
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=.5))+

ggplot(fix_libr_dt, aes(preparation, fill=count_lvl))+geom_bar(position="fill")+
  theme_classic()+facet_wrap(~library)+
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=.5))+
  scale_fill_manual(values = c("red4", "red1", "orange", "yellow",
                               "palegreen1","palegreen3", "darkgreen"))
ggsave("/DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s6_library_runs/are_the_missing_TAAs_the_same/general_stats.pdf",
       height=5, width=10)

lost=fix_libr_dt[counts==0,][,.N, by=c("library", "preparation")]
ggplot(lost, aes(preparation, fill=library, y=N))+geom_col(position = "dodge")+
  theme_classic()+theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=.5))

# 2/ Are the same ones problematic?
# ==================
lost_dt=split(fix_libr_dt[counts==0,], by=c("library", "preparation"))

# set for all venns
venn_theme <- theme(
  plot.margin = margin( r = 20, l = 20, unit = "mm")
)

# # big
lost_dt_big=lost_dt[grepl("big", names(lost_dt))]
names(lost_dt_big)=stringr::str_split_i(names(lost_dt_big),"\\.",2)

lost_dt_big_barcodes=lapply(lost_dt_big, function(x){x$barcode})
p1=ggVennDiagram::ggVennDiagram(lost_dt_big_barcodes,
                             label = "count", label_color = "white")+
  ggtitle("BIG_TAA: Lost barcode overlap")+
  venn_theme

lost_dt_big_genes=lapply(lost_dt_big, function(x){x$gene_name})
p2=ggVennDiagram::ggVennDiagram(lost_dt_big_genes,
                             label = "count", label_color = "white")+
  ggtitle("BIG_TAA: Lost gene overlap")+
  venn_theme

# # small
lost_dt_small=lost_dt[grepl("small", names(lost_dt))]
names(lost_dt_small)=stringr::str_split_i(names(lost_dt_small),"\\.",2)

lost_dt_small_barcodes=lapply(lost_dt_small, function(x){x$barcode})
p11=ggVennDiagram::ggVennDiagram(lost_dt_small_barcodes,
                                label = "count", label_color = "white")+
  ggtitle("small_TAA: Lost barcode overlap")+
  venn_theme

lost_dt_small_genes=lapply(lost_dt_small, function(x){x$gene_name})
p22=ggVennDiagram::ggVennDiagram(lost_dt_small_genes,
                                label = "count", label_color = "white")+
  ggtitle("small_TAA: Lost gene overlap")+
  venn_theme

pdf("/DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s6_library_runs/are_the_missing_TAAs_the_same/lost.overlaps.pdf", width = 12,height = 8)
cowplot::plot_grid(plotlist = c(p1,p2,p11,p22), align="hv", nrow = 2, greedy = T)
dev.off()

# across all
# -------------------

heat_libr_dt=dcast(fix_libr_dt, barcode + gene_name~ preparation+ library, value.var = "counts")
any0=apply(heat_libr_dt[, 3:ncol(heat_libr_dt), with=F],1,function(x){any(x==0, na.rm=T)})
dt=heat_libr_dt[any0,]

library(ComplexHeatmap)
library(circlize)
# BIG
# --- 2) build numeric matrix ---
cols <- names(heat_libr_dt)[grepl("big",names(heat_libr_dt), ignore.case = T )]
mat <- as.matrix(dt[, ..cols])
rownames(mat) <- paste(dt$barcode, dt$gene_name, sep = "_")

# --- 5) heatmap ---
pdf("/DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s6_library_runs/are_the_missing_TAAs_the_same/big_taa_heatmap_of_lost_barcodes.pdf",
    height=12, width = 4)
Heatmap(
  log10(na.omit(mat)+1),
  name = "counts",
  na_col = "white",
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  row_names_gp    = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  show_row_names = TRUE,
  column_names_rot = 90,
  
  row_title = "Barcodes with ≥1 zero",
  column_title = "Preparation × Library"
)
dev.off()


# SMALL
# --- 2) build numeric matrix ---
cols <- names(heat_libr_dt)[grepl("small",names(heat_libr_dt), ignore.case = T )]
mat <- as.matrix(dt[, ..cols])
rownames(mat) <- paste(dt$barcode, dt$gene_name, sep = "_")

# --- 5) heatmap ---

pdf("/DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s6_library_runs/are_the_missing_TAAs_the_same/small_taa_heatmap_of_lost_barcodes.pdf",
    height=12, width = 4)
Heatmap(
  log10(na.omit(mat)+1),
  name = "counts",
  na_col = "white",
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  column_names_rot = 90,
  
  row_title = "Barcodes with ≥1 zero",
  column_title = "Preparation × Library",
  row_names_gp    = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
)
dev.off()

# BOTH
# --- 2) build numeric matrix ---
cols <- names(heat_libr_dt)[3:ncol(heat_libr_dt)]
mat <- as.matrix(dt[, ..cols])
rownames(mat) <- paste(dt$barcode, dt$gene_name, sep = "_")

mat_bin <- matrix(NA_integer_, nrow = nrow(mat), ncol = ncol(mat),
                  dimnames = dimnames(mat))

mat_bin[!is.na(mat) & mat == 0]            <- 1L
mat_bin[!is.na(mat) & mat >= 1 & mat < 100] <- 2L
mat_bin[!is.na(mat) & mat >= 100 & mat < 1000] <- 3L
mat_bin[!is.na(mat) & mat >= 1000]           <- 4L
# NA in mat stays NA in mat_bin

col_fun <- c("0" = "grey70",    # NA in original
             "1" = "#2166ac",   # zero counts
             "2" = "lightyellow2",   # low counts [1,100)
             "3" = "#f4a582",
             "4"= "#d6604d")   # high counts ≥100

pdf("/DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s6_library_runs/are_the_missing_TAAs_the_same/both_taa_heatmap_of_lost_barcodes.pdf",
    height=12, width = 4)
Heatmap(
  mat_bin,
  name = "counts",
  na_col = "grey70",
  col = col_fun,
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  column_names_rot = 90,
  row_title = "Barcodes with ≥1 zero",
  column_title = "Preparation × Library",
  row_names_gp    = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    at = c(1, 2, 3,4),
    labels = c("zero", "[1,100)","[100,1000)", "≥1000")
  )
)
dev.off()
