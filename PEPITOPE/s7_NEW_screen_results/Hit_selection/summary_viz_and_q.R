library(data.table)
library(ggplot2)
library(ggrepel)
theme_set(theme_classic())
set.seed(150799)

output_dir <- "PEPITOPE/s7_NEW_screen_results/Hit_selection/summary_output/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

screen_1  <- fread("PEPITOPE/s7_NEW_screen_results/output_screen_1/updated_screen_results/screen1_results.txt")
screen_2a <- fread("PEPITOPE/s7_NEW_screen_results/output_screen_2a/updated_screen_results/screen2a_results.txt")

plot_summary=function(results){
  results_f=results[,.(gene_id, ref_key, mutation_profile, comparison, barcode, baseMean,
                     log2FoldChange, lfcSE, stat, pvalue, padj)]
  results_f=split(results_f, f=results_f$comparison)
  results_f=lapply(results_f, function(x){x[, comparison:=NULL]})
  comp_db_sg_ut=merge(results_f$`Cluster vs UT`, results_f$`Single vs UT`,
                      by=c("gene_id","barcode", "ref_key","mutation_profile"), 
                      suffixes=c("_db", "_sg"))
  comp_db_sg_ut[, mutation_profile:= ifelse(grepl("CMV|EBV|HPV|FLU|SARS", gene_id),
                                            "Viral", mutation_profile)]
  
  ggplot(comp_db_sg_ut[baseMean_db>100 | baseMean_sg>100], aes(stat_sg, stat_db))+
    geom_point(aes(col=mutation_profile), size = 3, alpha=1)+
    ggrepel::geom_text_repel(data=comp_db_sg_ut[mutation_profile!="ref" &
                                                (baseMean_db>100 | baseMean_sg>100) & 
                                                  (padj_db<0.05 | padj_sg<(0.05)) &
                                                  (log2FoldChange_db<(-0.5) | log2FoldChange_sg<(-0.5))],
                             aes(label=gene_id), size=4, max.overlaps = 25)+
    scale_color_manual(values=c(alt="red3", ref="grey", "TAA"="gold", "Viral"="steelblue2"))+
    labs(y="Stat: Cluster TCRs vs UT", x="Stat: Singlet TCRs vs UT")+
    theme(legend.position = "none", axis.title = element_text(size=14))+
    geom_hline(yintercept = 0, lty="dashed", col="grey")+
    geom_vline(xintercept = 0, lty="dashed", col="grey")
}

plot_summary(screen_1)
plot_summary(screen_2a)
