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

# we had previously sequenced the input libraries, albeit separately
# we want to see if there were any issues dragged from those
# eg barcodes that were lost from the beginning

# previous run -- TAA libraries
taa_lib_input_sample1=fread("PEPITOPE/s6_library_runs/second_run_may_26/output_sample_1/report_source_data_pat_lib_plus_bigTAA.txt")
taa_lib_input_sample2=fread("PEPITOPE/s6_library_runs/second_run_may_26/output_sample_2/report_source_data_pat_lib_plus_bigTAA.txt")

# previous run -- All libraries
all_lib_input_sample1=fread("PEPITOPE/s6_library_runs/first_run_march_26/output/report_source_data_pat_lib_plus_bigTAA.txt")
all_lib_input_sample1=all_lib_input_sample1[expected=="Expected"]

# current stats
new_screen_1=fread("PEPITOPE/s7_screen_results/output_screen_1/manual_qc_rows.txt")
new_screen_2a=fread("PEPITOPE/s7_screen_results/output_screen_2a/manual_qc_rows.txt")


# 1/ Compare screen counts with counts in the original libraries
new_screen_2a[, unmatched:= NULL]
new_screen_2a[, mean_guide:= rowMeans(.SD), by=c("guide", "gene", "guide_type")]

corr_dt=merge(all_lib_input_sample1,new_screen_2a, by="guide")
p1=ggplot(corr_dt, aes(counts+1, mean_guide+1))+
  geom_point(aes(col=patient_id), alpha=.3, size=1, stroke=.5)+
  ggtitle("Guide counts (+1)")+
  labs(y="Screen 2a avg. guide counts", x="Library 1st test (all barcodes)")+
  scale_x_log10()+scale_y_log10()

new_screen_1[, unmatched:= NULL]
new_screen_1[, mean_guide:= rowMeans(.SD), by=c("guide", "gene", "guide_type")]

corr_dt_sc1=merge(all_lib_input_sample1,new_screen_1, by="guide")
p2=ggplot(corr_dt_sc1, aes(counts+1, mean_guide+1))+
  geom_point(aes(col=patient_id), alpha=.3, size=1, stroke=.5)+
  ggtitle("Guide counts (+1)")+
  labs(y="Screen 1 avg. guide counts", x="Library 1st test (all barcodes)")+
  scale_x_log10()+scale_y_log10()

# 2/ Compare with the big TAA second library

corr_dt_2_sc2a=merge(taa_lib_input_sample2,new_screen_2a, by="guide")
p3=ggplot(corr_dt_2_sc2a, aes(counts+1, mean_guide+1))+
  geom_point(aes(col=patient_id), alpha=.3, size=1, stroke=.5)+
  ggtitle("Guide counts (+1)")+
  labs(y="Screen 2a avg. guide counts", x="Library 2st test (only TAA)")+
  scale_x_log10()+scale_y_log10()

corr_dt_2_sc1=merge(taa_lib_input_sample2,new_screen_1, by="guide")
p4=ggplot(corr_dt_2_sc1, aes(counts+1, mean_guide+1))+
  geom_point(aes(col=patient_id), alpha=.3, size=1, stroke=.5)+
  ggtitle("Guide counts (+1)")+
  labs(y="Screen 1 avg. guide counts", x="Library 2st test (only TAA)")+
  scale_x_log10()+scale_y_log10()


cowplot::plot_grid(plotlist = list(p1,p2,p3,p4))
ggsave("PEPITOPE/s7_screen_results/output_s3/input_library_to_measured_counts.pdf",
       width = 7, height=5, create.dir = T)
