library(pepitope)
library(data.table)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
set.seed(150799)

# 1/ Sample demultiplexing 
# ================================================

sample_sheet=data.table(
  sample_id=c("P1", "P2", "P3", "P4", "P5", "TAA_big", "TAA_small"),
  patient=c("P1", "P2", "P3", "P4", "P5", "TAA_big", "TAA_small"),
  rep=rep(1, 7),
  origin=c("neo", "neo", "neo", "neo", "neo", "TAA", "TAA"),
  barcode=c("AAGACCA", "CCAGTGT", "TGAGTCC","CAAGATG","AACCGAC","AGAATCG","AGACCGT"))

# # needed to install fqtk in a conda environment
Sys.setenv(PATH = paste(
  "/home/j.simon/miniconda3/envs/fqtk_env/bin",
  Sys.getenv("PATH"),
  sep = ":"))

Sys.which("fqtk")
temp_dir=demux_fq("PEPITOPE/s6_library_runs/first_run_march_26/8613_1_BQC001_GAATTGCGAA-CCGGCACTAT_S32_R1_001.fastq.gz",
                  samples = sample_sheet, read_structures = "7B+T")

if(!dir.exists("PEPITOPE/s6_library_runs/first_run_march_26/output/")){
  dir.create("PEPITOPE/s6_library_runs/first_run_march_26/output/")}

# # visualise results
demux_metrics=fread(paste0(temp_dir, "/demux-metrics.txt"))
dmm=melt(demux_metrics, id.vars=c("sample_id", "barcode"))
ggplot(dmm, aes(sample_id, value))+
  facet_wrap(~variable, scale="free_y", nrow=1)+
  geom_col(fill="coral3")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
ggsave("PEPITOPE/s6_library_runs/first_run_march_26/output/sanity_metrics.pdf",
       height=4, width = 10)


# 2/ CREATE FUNCTIONS INPUT LIST:
# ============================================

# do we need reverse complements?
reverse_complement=TRUE
# -----

# # Grab patient barcodes
barcode_dt=read_xlsx("PEPITOPE/s5_create_libraries/resources_zm/minigene_barcode_list_all.xlsx")
early_barcodes=barcode_dt$Barcode

# # Grab TAA barcodes
taa_big=fread("PEPITOPE/s6_library_runs/shared_by_JV/large_TAA.csv")
taa_small=fread("PEPITOPE/s6_library_runs/shared_by_JV/small_TAA_fix.csv")
# # reconfigs
setnames(taa_big, "BARCODE", "barcode")
taa_big=taa_big[!is.na(LENGTH)]
setnames(taa_small, "Oligo barcode", "barcode")
taa_small=taa_small[!duplicated(taa_small)]
# exact repetition of the same row? 4764 EBV_EBNA1   ACTGTCA

# # all barcodes:
ext_valid_barcodes=union(early_barcodes, union(taa_big$barcode, taa_small$barcode))

# #  Create library .tsv
tsv = tibble(name = as.character(ext_valid_barcodes)) |>
mutate(barcode=name, gene=name)

if (reverse_complement)
  tsv$barcode = as.character(reverseComplement(DNAStringSet(tsv$barcode)))
lpath = file.path(temp_dir, "lib.tsv")
utils::write.table(tsv, file=lpath, sep="\t", row.names=FALSE, quote=FALSE)

# # relocate results
if(!dir.exists("PEPITOPE/s6_library_runs/first_run_march_26/tmp_files/")){
  dir.create("PEPITOPE/s6_library_runs/first_run_march_26/tmp_files/")}
R.utils::copyDirectory( from=temp_dir, 
                        to="PEPITOPE/s6_library_runs/first_run_march_26/tmp_files/", 
                        recursive=T)

# 2.1 Reformat the minigene libraries for downstream plots
# ============================================
all_constructs=fread("PEPITOPE/s5_create_libraries/new_rna_library_all_patients_duplo.txt")
setnames(all_constructs, "MINIGENE", "tiled")
setnames(all_constructs, "BARCODE", "barcode")
pat_dt=all_constructs[, .(barcode, replicate, patient_id, gene_name, mut_id, pep_id, pep_type)]


# *rename big list
setnames(taa_big, "GENE NAME", "gene_name")
setnames(taa_big, "MINIGENE SEQ", "tiled")
taa_big[, pep_type:= "TAA"]
taa_big[, mut_id:= "TAA"]
taa_big[, pep_id:= "TAA"]
taa_big[, sample_id:= "TAA_big"]
taa_big[, patient_id:= "TAA_big"]

# *rename small list
setnames(taa_small, "GENE NAME", "gene_name")
setnames(taa_small, "MINIGENE SEQ", "tiled")
taa_small[, pep_type:= "TAA"]
taa_small[, mut_id:= "TAA"]
taa_small[, pep_id:= "TAA"]
taa_small[, sample_id:= "TAA_small"]
taa_small[, patient_id:= "TAA_small"]

# create patients + TAA big
tot_barcode_anno_big=rbind(
  taa_big[, .(barcode, replicate="1", patient_id, gene_name, mut_id, pep_id, pep_type)], 
  pat_dt)

unused_or_taa_small=setdiff(ext_valid_barcodes, tot_barcode_anno_big$barcode)
filler_big=data.table(barcode=unused_or_taa_small, replicate="NS", 
                      patient_id="NS", gene_name="NS", mut_id="NS", 
                      pep_id="NS", pep_type="NS")

tot_barcode_anno_big_filled=rbind(tot_barcode_anno_big, filler_big)

# create patients + TAA big
tot_barcode_anno_small=rbind(
  taa_small[, .(barcode, replicate="1", patient_id, gene_name, mut_id, pep_id, pep_type)], 
  pat_dt)

unused_or_taa_big=setdiff(ext_valid_barcodes, tot_barcode_anno_small$barcode)
filler_small=data.table(barcode=unused_or_taa_big, replicate="NS", 
                      patient_id="NS", gene_name="NS", mut_id="NS", 
                      pep_id="NS", pep_type="NS")

tot_barcode_anno_small_filled=rbind(tot_barcode_anno_small, filler_small)


# 3/ RUN BASH FUNCTION manual_guide_counter_run.sh
# ============================================
system('bash /DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s6_library_runs/first_run_march_26/manual_guide_counter_run.sh')

# 4/ GET RESULTS FROM GUIDE COUNTER
# ============================================
new_tmp="PEPITOPE/s6_library_runs/first_run_march_26/tmp_files/"
manual_output=list.files(new_tmp, pattern = "manual", full.names = T)

stats=fread(manual_output[which(grepl("stat",manual_output))])
counts=fread(manual_output[which(grepl("manual_barcodes.counts.txt",manual_output))])
rows=fread(manual_output[which(grepl("extended",manual_output))])

stats = stats |> mutate(sample_id = sub("\\.R1$", "", label))
colnames(counts) = sub("\\.R1$", "", colnames(counts))
colnames(rows) = sub("\\.R1$", "", colnames(rows))


# 5/ BUILD FINAL REPORT
# ============================================

# # splitting TAA libraries --> BIG
row_meta_big=merge(rows[, barcode:= guide], tot_barcode_anno_big_filled, by="barcode", all=T)
fwrite(row_meta_big, "PEPITOPE/s6_library_runs/first_run_march_26/output/res_pat_lib_plut_bigTAA.txt")
m_big=melt(row_meta_big, 
           measure.vars = c("P1","P2","P3","P4","P5","TAA_big"), 
           variable.name = "detected_in", value.name = "counts")

# # splitting TAA libraries --> SMALL
row_meta_small=merge(rows, tot_barcode_anno_small_filled, by="barcode", all=T)
fwrite(row_meta_small, "PEPITOPE/s6_library_runs/first_run_march_26/output/res_pat_lib_plut_smallTAA.txt")
m_small=melt(row_meta_small, 
           measure.vars = c("P1","P2","P3","P4","P5","TAA_small"), 
           variable.name = "detected_in", value.name = "counts")


plot_results=function(melt_dt, filename){
  
  melt_dt[, expected:= factor(ifelse(patient_id==detected_in, "Expected", "Contamination"))]
  melt_dt[, count_lvl:= cut(counts, breaks = c(0,1,10,100,300,500,1000,Inf), include.lowest = T, ordered_result = T)]
  melt_dt[, rank:= rank(counts, ties.method = "random"), by=c("detected_in", "expected")]
  melt_dt[, rank_norm:= rank/.N, by=c("detected_in", "expected")]
  totals_dt <- melt_dt[, .(total = sum(counts)), by = detected_in]
  fwrite(melt_dt, paste0("PEPITOPE/s6_library_runs/first_run_march_26/output/report_source_data_",filename,".txt"))
  
  melt_dt[, counts_non_0:= counts>0]
  test_recovery=melt_dt[, .N, by=c("expected", "detected_in", "counts_non_0")]
  test_recovery=test_recovery[expected=="Expected"]
  
  p1=ggplot(melt_dt, aes(detected_in, y=counts, fill=expected))+
    geom_col(position = "stack")+theme_classic()+ggtitle("Total counts")+
    geom_text(
      data = totals_dt,
      aes(x = detected_in, y = total, label = round(total, 1)),
      inherit.aes = FALSE,
      vjust = -0.5,
      size = 3
    )
  p2=ggplot(melt_dt, aes(detected_in, y=counts, fill=expected))+geom_col(position = "fill")+theme_classic()+ggtitle("Total counts distribution")
  p3=ggplot(test_recovery, aes(detected_in, y=N, fill=counts_non_0))+geom_col()+
    geom_text(aes(label = N),
              position = position_stack(vjust = 0.5),size = 3.5) +theme_classic() +
    ggtitle("Minigenes in library (detected?)")+
    scale_fill_manual(values=c("palegreen3", "coral1"))
  p4=ggplot(melt_dt, aes(detected_in, fill = expected, y=counts+1))+geom_boxplot()+
    scale_y_log10()+theme_classic()+ggtitle("Median detection rate") +
    stat_summary(
      fun = function(y) median(y),        # median of (counts+1)
      geom = "text",
      aes(label = round(after_stat(y) - 1, 1)),  # subtract 1 to show original scale
      position = position_dodge(width = 0.75),
      vjust = -1.5,
      size = 3
    ) 
  p5=ggplot(melt_dt, aes(expected,fill=count_lvl))+geom_bar(position='fill')+
    theme_classic()+facet_wrap(~detected_in, nrow=1)+
    theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=.5))+
    scale_fill_manual(values = c("red4", "red1", "orange", "yellow","palegreen1","palegreen3", "darkgreen"))+
    ggtitle("Detection rate groups distribution")
  p6=ggplot(melt_dt, aes(rank_norm, counts+1, col=detected_in))+geom_point(aes(shape=expected), alpha=.4)+
    geom_hline(yintercept = c(10,100,1000), lty="dashed", col="grey")+
    geom_vline(xintercept = c(0.05), lty="dashed", col="grey")+
    theme_classic()+
    scale_y_continuous(breaks=c(0,10,100,1000,10000), transform = "log10")+
    scale_x_continuous(breaks=c(0,0.05,0.5,1))+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
    ggtitle("Per barcode/library: Detection rate distribution")


cowplot::plot_grid(plotlist = list(p1,p2,p3,p4,p5,p6), align = "hv", ncol = 2, axis = "lrtb")
ggsave(paste0("PEPITOPE/s6_library_runs/first_run_march_26/output/report_", filename, ".pdf"),
    width=12, height=18)
}

plot_results(m_big, "pat_lib_plus_bigTAA")
plot_results(m_small, "pat_lib_plus_smallTAA")
