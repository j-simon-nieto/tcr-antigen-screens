library(data.table)
library(ggplot2)

set_theme(theme_classic())
set.seed(150799)
output_dir="PEPITOPE/s7_screen_results/answering_JV_observations/unmatched/"
if(!dir.exists(output_dir)){dir.create(output_dir)}

# Load data
counts_screen_1=fread("PEPITOPE/s7_screen_results/output_screen_1/manual_qc_rows.txt")
vars=names(counts_screen_1)[4:ncol(counts_screen_1)]
counts_screen_1[, tot_barcode_reads:= rowSums(.SD), .SDcols = vars, by=c("guide", "gene")]
counts_screen_1[, unmatched_fraction:= 100*unmatched/tot_barcode_reads]
counts_screen_1=counts_screen_1[tot_barcode_reads>0]

# Get reference
tot_barcode=fread("PEPITOPE/s7_screen_results/answering_JV_observations/missing_barcodes/tot_barcode.txt") # from prev script

# most of the unmatched reads come from barcodes with low reads anyway (<100)
ggplot(counts_screen_1, aes(tot_barcode_reads, unmatched_fraction))+
  geom_point(alpha=.2, size=.2, position = position_jitter(width = .1, height = 2))+
  scale_x_log10(breaks=c(1,10,50,100,500,1000,5000,10000,50000, 100000))+
  geom_point(size=.5, col="red")
ggsave(paste0(output_dir, "only_low_count_barcodes_have_preponderant_unamtched_populations.pdf"), width=5, height = 4)

# are there differences?
counts_screen_1[, guide_type:= stringr::str_split_i(gene, "_", -2)]
counts_screen_1[, unmatched_percent_bins:= cut(unmatched_fraction, 
                                               breaks=c(0,2,100), 
                                               include.lowest=T)]

ggplot(counts_screen_1, aes(unmatched_percent_bins, fill=guide_type))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = c(alt = "red3", ref = "grey20", TAA = "gold3"))+
  xlab("% of unmatched reads in barcode")

# contamination?
pat="P1"
counts_screen_1[, good_cont:= ifelse(guide %in% 
                                       tot_barcode[patient_id %in% c("TAA_big", pat)]$barcode,
                "good", "contamination")]

ggplot(counts_screen_1, aes(unmatched_percent_bins, fill=guide_type))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = c(alt = "red3", ref = "grey20", TAA = "gold3"))+
  xlab("% of unmatched reads in barcode")+
  facet_wrap(~good_cont)+
  ggplot(counts_screen_1, aes(unmatched_percent_bins, fill=guide_type))+
  geom_bar(position = "stack")+
  scale_fill_manual(values = c(alt = "red3", ref = "grey20", TAA = "gold3"))+
  xlab("% of unmatched reads in barcode")+
  facet_wrap(~good_cont)
ggsave(paste0(output_dir, "unmatched_dsitribution_for_P1_is_uniform.pdf"), width=8, height = 4)


ggplot(counts_screen_1, aes(tot_barcode_reads, unmatched_fraction))+
  geom_point(alpha=.2, size=.2, position = position_jitter(width = .1, height = 2))+
  scale_x_log10(breaks=c(1,10,50,100,500,1000,5000,10000,50000, 100000))+
  geom_point(size=.5, col="red")+facet_wrap(~good_cont, nrow=2)+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
ggsave(paste0(output_dir, "remove_contaminations_only_low_count_barcodes_have_preponderant_unamtched_populations.pdf"), width=5, height = 6)
