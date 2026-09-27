library(pepitope)
library(data.table)
library(ggplot2)
library(dplyr)
library(SummarizedExperiment)
library(readxl)
library(Biostrings)
set.seed(150799)

# Patient 5 (SJ32), standard screen 

output_dir="PEPITOPE/s7_NEW_screen_results/output_screen_2a/"
if(!dir.exists(output_dir)){
  dir.create(output_dir)}

exp_tmp_dir="PEPITOPE/s7_NEW_screen_results/tmp_files_screen_2a/"

# Grabbing TAA big library for barcodes

# 0/ Load files
# ================================================
sample_info=fread("PEPITOPE/s7_screen_results/input/metadata/screen_2_SJ32_primer_pcr_data.csv")[patient_id=="SJ32"]
sample_info=sample_info[!is.na(`Total volume`)]


# 1/ Sample demultiplexing 
# ================================================

sample_sheet=data.table(
  sample_id=sample_info$sample_contains,
  patient=sample_info$patient_id,
  rep=stringr::str_split_i(sample_info$sample_contains, " - ", 2),
  origin=stringr::str_split_i(sample_info$sample_contains, " - ", 1), 
  primer_1=sample_info$`Primer Fw`,
  barcode=sample_info$fw_primer_seq,
  barcode_rev_complement=sample_info$rev_comp_primer_seq,
  primer2=sample_info$`Primer Rv`,
  rev_complement=sample_info$rev_comp_primer_seq
)

# # needed to install fqtk in a conda environment
Sys.setenv(PATH = paste(
  "/home/j.simon/miniconda3/envs/fqtk_env/bin",
  Sys.getenv("PATH"),
  sep = ":"))

Sys.which("fqtk")
# # NOTE: screen2a only
temp_dir=demux_fq("PEPITOPE/s7_screen_results/input/fastq_files/8793_2_BSCREEN002-a_TCCTGGAACA-TACTGGTTGT_S2_R1_001.fastq.gz",
                  samples = sample_sheet, read_structures = "7B+T")

# # visualise results
demux_metrics=fread(paste0(temp_dir, "/demux-metrics.txt"))
dmm=melt(demux_metrics, id.vars=c("sample_id", "barcode"))

palette_sample_id=c("1D3 - A"="gold1","1D3 - B"="gold3","1D3 - C"="gold4",
                    "B only - A"="grey80","B only - B"="grey70","B only - C"="grey60",
                    "Cluster T - A"="red1","Cluster T - B"="red3","Cluster T - C"="red4",
                    "Single T - A"="steelblue1","Single T - B"="steelblue3","Single T - C"="steelblue4", 
                    "UT - A"="pink1","UT - B"="pink3","UT - C"="pink4",
                    "unmatched"="black")

ggplot(dmm, aes(sample_id, value))+
  facet_wrap(~variable, scale="free_y", nrow=1)+
  geom_col(aes(fill=sample_id))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
  scale_fill_manual(values=palette_sample_id)
ggsave(paste0(output_dir,"/sanity_metrics.pdf"),
       height=4, width = 10)


# 2/ CREATE FUNCTIONS INPUT LIST:
# ============================================

# do we need reverse complements?
reverse_complement=TRUE
# -----

# # Grab TAA barcodes
taa_big=fread("PEPITOPE/s6_library_runs/shared_by_JV/large_TAA.csv")

# # reconfigs
setnames(taa_big, "BARCODE", "barcode")
taa_big=taa_big[!is.na(LENGTH)]

# *rename big list
setnames(taa_big, "GENE NAME 2", "gene_name")
setnames(taa_big, "MINIGENE SEQ", "tiled")
taa_big[, pep_type:= "TAA"]
taa_big[, mut_id:= "TAA"]
taa_big[, pep_id:= "TAA"]
taa_big[, sample_id:= "TAA_big"]
taa_big[, patient_id:= "TAA_big"]

# 2.1 Reformat the minigene libraries for downstream plots
# ============================================
all_constructs=fread("PEPITOPE/s5_create_libraries/new_rna_library_all_patients_duplo.txt")
setnames(all_constructs, "MINIGENE", "tiled")
setnames(all_constructs, "BARCODE", "barcode")
pat_dt=all_constructs[, .(barcode, replicate, patient_id, gene_name, mut_id, pep_id, pep_type)]

# create patients + TAA big
tot_barcode_anno_big=rbind(
  taa_big[, .(barcode, replicate="1", patient_id, gene_name, mut_id, pep_id, pep_type)], 
  pat_dt)

tot_barcode_anno_big_filled=tot_barcode_anno_big

tot_barcode_anno_big_filled[, gene_id:= ifelse(pep_id=="TAA", stringr::str_split_i(gene_name,"-",1),
                                                 ifelse(gene_name=="", gsub("\\.fs", "",mut_id),gene_name))]

# # all barcodes:
ext_valid_barcodes=c(tot_barcode_anno_big_filled$barcode)

# #  Create library .tsv --> supervise
tsv = tibble(name= as.character(tot_barcode_anno_big_filled$barcode), #otherwise it breaks
             barcode = as.character(tot_barcode_anno_big_filled$barcode),
             gene=paste(tot_barcode_anno_big_filled$gene_name,
                        tot_barcode_anno_big_filled$pep_id,
                        tot_barcode_anno_big_filled$pep_type,
                        tot_barcode_anno_big_filled$replicate,sep = "_"))

if(reverse_complement)
  tsv$barcode = as.character(reverseComplement(DNAStringSet(tsv$barcode)))
lpath = file.path(temp_dir, "lib.tsv")
utils::write.table(tsv, file=lpath, sep="\t", row.names=FALSE, quote=FALSE)

# # relocate results
if(!dir.exists(exp_tmp_dir)){dir.create(exp_tmp_dir)}
R.utils::copyDirectory( from=temp_dir, 
                        to=exp_tmp_dir, 
                        recursive=T)


# 3/ RUN BASH FUNCTION manual_guide_counter_run.sh
# ============================================

# # change accordingly
run=paste0("guide-counter count ",
           "--offset-min-fraction 0.2 " ,
           "--input /DATA/j.simon/JV/ANTIGEN_calling/",exp_tmp_dir,"*.R1.fq.gz ", # <-- needs to grab the correct ones
           "--library /DATA/j.simon/JV/ANTIGEN_calling/",exp_tmp_dir,"lib.tsv ",
           "-x ",
           "-o /DATA/j.simon/JV/ANTIGEN_calling/",exp_tmp_dir,"manual_barcodes")
system(run)

# 4/ GET RESULTS FROM GUIDE COUNTER
# ============================================
manual_output=list.files(exp_tmp_dir, pattern = "manual", full.names = T)

stats=fread(manual_output[which(grepl("stat",manual_output))])
counts=fread(manual_output[which(grepl("manual_barcodes.counts.txt",manual_output))])
rows=fread(manual_output[which(grepl("extended",manual_output))])

stats = stats |> mutate(sample_id = sub("\\.R1$", "", label))
colnames(counts) = sub("\\.R1$", "", colnames(counts))
colnames(rows) = sub("\\.R1$", "", colnames(rows))

fwrite(stats, paste0(output_dir, "manual_qc_stats.txt"))
fwrite(counts, paste0(output_dir, "manual_qc_counts.txt"))
fwrite(rows, paste0(output_dir, "manual_qc_rows.txt"))


# 5/ BUILD FINAL REPORT
# ============================================

# # investigating only the big TAA
row_meta_big=merge(rows[, barcode:= guide], tot_barcode_anno_big_filled, by="barcode", all=T)
fwrite(row_meta_big, paste0(output_dir,"res_pat_lib_plut_bigTAA.txt"))
m_big=melt(row_meta_big, 
           measure.vars = c(sample_info$sample_contains, "unmatched"), # requires attention
           variable.name = "detected_in", value.name = "counts")
m_big_f=m_big[  !is.na(counts),]


plot_results=function(melt_dt, filename){
  
  melt_dt[, expected:= factor(ifelse(patient_id %in% c("P5","TAA_big")| pep_type %in% c("TAA"), "Expected", "Contamination"))] # contamination are barcodes from other patients
  melt_dt[, count_lvl:= cut(counts, breaks = c(0,1,10,100,300,500,1000,Inf), include.lowest = T, ordered_result = T)]
  melt_dt[, rank:= rank(counts, ties.method = "random"), by=c("detected_in", "expected")]
  melt_dt[, rank_norm:= rank/.N, by=c("detected_in", "expected")]
  totals_dt <- melt_dt[, .(total = sum(counts)), by = detected_in]
  fwrite(melt_dt, paste0(output_dir,"report_source_data_",filename,".txt"))
  
  melt_dt[, counts_non_0:= counts>0]
  test_recovery=melt_dt[, .N, by=c("expected", "detected_in", "counts_non_0")]
  test_recovery=test_recovery[expected=="Expected"]
  
  p1=ggplot(melt_dt, aes(detected_in, y=counts, fill=expected))+
    geom_col(position = "stack")+theme_classic()+ggtitle("Total counts")+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
    geom_text(
      data = totals_dt,
      aes(x = detected_in, y = total, label = round(total, 1)),
      inherit.aes = FALSE,
      vjust = -0.5,
      size = 3
    )
    
  p2=ggplot(melt_dt, aes(detected_in, y=counts, fill=expected))+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
      geom_col(position = "fill")+theme_classic()+ggtitle("Total counts distribution")
    
  p3=ggplot(test_recovery, aes(detected_in, y=N, fill=counts_non_0))+geom_col()+
    geom_text(aes(label = N),
              position = position_stack(vjust = 0.5),size = 3.5) +theme_classic() +
    ggtitle("Minigenes in library (detected?)")+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
    scale_fill_manual(values=c("palegreen3", "coral1"))
  
  p4=ggplot(melt_dt, aes(detected_in, fill = expected, y=counts+1))+geom_boxplot()+
    scale_y_log10()+theme_classic()+ggtitle("Median detection rate") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
    stat_summary(
      fun = function(y) median(y),        # median of (counts+1)
      geom = "text",
      aes(label = round(after_stat(y) - 1, 1)),  # subtract 1 to show original scale
      position = position_dodge(width = 0.75),
      vjust = -1.5,
      size = 3
    ) 
  p5=ggplot(melt_dt, aes(detected_in,fill=count_lvl))+geom_bar(position='fill')+
    theme_classic()+facet_wrap(~expected, nrow=1)+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
    scale_fill_manual(values = c("red4", "red1", "orange", "yellow","palegreen1","palegreen3", "darkgreen"))+
    ggtitle("Detection rate groups distribution")
  
  p6=ggplot(melt_dt, aes(rank_norm, counts+1, col=detected_in))+geom_point(aes(shape=expected), alpha=.4)+
    geom_hline(yintercept = c(10,100,1000), lty="dashed", col="grey")+
    geom_vline(xintercept = c(0.05), lty="dashed", col="grey")+
    theme_classic()+
    scale_y_continuous(breaks=c(0,10,100,1000,10000), transform = "log10")+
    scale_x_continuous(breaks=c(0,0.05,0.5,1))+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
    scale_color_manual(values=palette_sample_id)+
    ggtitle("Per barcode/library: Detection rate distribution")
  
  
  cowplot::plot_grid(plotlist = list(p1,p2,p3,p4,p5,p6), align = "hv", ncol = 2, axis = "lrtb")
  ggsave(paste0(output_dir,"/report_", filename, ".pdf"),
         width=12, height=18)
}

plot_results(m_big_f, "pat_lib_plus_bigTAA")
