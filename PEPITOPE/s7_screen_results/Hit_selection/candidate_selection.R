# GOALS:
# -------------------------
# select which hits from each screen to move forward with

# Main criteria:
# -------------------------
# significant difference in reads
# TAA: strong outlier + average population
# Neoantigen: prioritise at least one strong outlier and no ref effects
# Viral (manual annotation): remove donor specificities (drop outs in UT v B)
# Prioritise 10-15 per patient (DB+SG)

library(data.table)
library(ggplot2)
library(ggrepel)
set_theme(theme_classic())
set.seed(150799)
output_dir="PEPITOPE/s7_screen_results/Hit_selection/output/"
if(!dir.exists(output_dir)){dir.create(output_dir)}

# Load data
screen_1=fread("PEPITOPE/s7_screen_results/output_screen_1/screen_results/screen1_results.txt")
screen_2a=fread("PEPITOPE/s7_screen_results/output_screen_2a/screen_results/screen2a_results.txt")

# Create function
hit_prioritisation=function(screen_dt, minBase=100,
                            COMP_MAIN,COMP_CTRL="UT vs B",
                            folder_output){
  # # filter to comparison of interest and reliable baseMean
  screen_dt=screen_dt[baseMean>minBase,]
  screen_dt=screen_dt[, if(.N > 1) .SD, by=c("gene_id","comparison")] # at least found twice per comp (allowing only one ref or one alt)
  setorder(screen_dt, stat)
  screen_dt[, min_stat:= min(stat), by=c("gene_id","mutation_profile", "comparison")] # best replicate
  screen_dt[, max_stat:= max(stat), by=c("gene_id","mutation_profile", "comparison")] # worse replicate
  
  screen_dt[, sd_x:=sd(stat), by="comparison"]
  screen_dt[, mean_x:=mean(stat), by="comparison"]
  screen_dt[, thr:=mean_x-sd_x ] # one standard deviation over mean
  screen_dt[, thr2:=mean_x-2*sd_x] # two sds over mean
  
  # main comparison:
  screen_db_ut=screen_dt[comparison==COMP_MAIN]
  
  # focus on TAAs
  sc1_taa=screen_db_ut[mutation_profile=='TAA']
  sc1_taa=sc1_taa[ stat <=thr & padj<0.05 & max_stat<0] # # # #THIS DECISION PROPAGATES: can either be thr & 0 or thr2 and thr
  fwrite(sc1_taa, paste0(folder_output, gsub(" ","_",COMP_MAIN), "__TAA_hits_stats.txt"))
  final_taa=sc1_taa[,.N, by=c("gene_id","min_stat")]
  
  # control: avoid donor derived biases (should already be implicit in cluster vs UT/1d3)
  screen_control=screen_dt[comparison==COMP_CTRL]
  sc1_taa_ctrl=screen_control[mutation_profile=='TAA']
  sc1_taa_ctrl_f=sc1_taa_ctrl[ stat <=thr & padj<0.05 & max_stat<0] # # # #THIS DECISION PROPAGATES: can either be thr & 0 or thr2 and thr
  fwrite(sc1_taa_ctrl_f, paste0(folder_output, gsub(" ","_",COMP_MAIN), "__TAA_control_stats.txt"))
  final_ctrl_taa=sc1_taa_ctrl_f[,.N, by=c("gene_id","min_stat")]
  
  compare_taa=merge(final_taa,unique(sc1_taa_ctrl[,.(gene_id, min_stat)]),
                    suffixes=c("", "_ctrl"), 
                    by="gene_id", all.x=T)
  
  compare_taa[, was_hit_in_ctrl:= gene_id %in% final_ctrl_taa$gene_id]
  fwrite(compare_taa, paste0(folder_output, gsub(" ","_",COMP_MAIN), "__TAA_hits.txt"))
  
  # focus on NeoAntigens
  screen_db_ut[, min_stat_ctrl:= min(stat[mutation_profile=='ref']), by="gene_id"]
  screen_db_ut[, mean_ref:= mean(stat[mutation_profile=='ref']), by="gene_id"]
  screen_db_ut[, mean_alt:= mean(stat[mutation_profile=='alt']), by="gene_id"]
  screen_db_ut[, alt_ref:= mean_alt-mean_ref<1] # arbitrary cutoff to avoid very marginal numbers
  
  sc1_alt=screen_db_ut[mutation_profile=='alt']
  sc1_alt=sc1_alt[ min_stat <=thr & padj<0.05 & alt_ref] # # # #THIS DECISION PROPAGATES: can either be thr2 or thr
  fwrite(sc1_alt, paste0(folder_output, gsub(" ","_",COMP_MAIN), "__neoantigen_hits_stats.txt"))
  final_alt=sc1_alt[,.N,  by=c("gene_id","min_stat","min_stat_ctrl")]
  fwrite(final_alt, paste0(folder_output, gsub(" ","_",COMP_MAIN), "__neoantigen_hits.txt"))
  
  # cross lists
  final_hits=rbindlist(list("TAA"=compare_taa,"ALT"=final_alt), idcol = "origin", use.names = T, fill = TRUE)
  setorder(final_hits, min_stat)
  final_hits[, rank:= frank(min_stat)]
  final_hits[, margin:= min_stat_ctrl-min_stat]
  fwrite(final_hits, paste0(folder_output, gsub(" ","_",COMP_MAIN), "__all_hits.txt"))
  
  return(final_hits)
}

create_list=function(screen_clusters, screen_singlets, MARGIN=2,
                     folder_output){
  screen_summary=merge(screen_clusters[,screen:="DB"],
                         screen_singlets[,screen:="SG"], 
                         by=c("origin", "gene_id", "was_hit_in_ctrl"),
                         suffixes=c("_DB", "_SG"), all = T)
  
  screen_summary[, min_stat_overall:=pmin(min_stat_DB,min_stat_SG, na.rm=T)]
  screen_summary[, margin:= pmax(margin_DB,margin_SG, na.rm=T)]
  setorder(screen_summary, min_stat_overall)
  fwrite(screen_summary, paste0(folder_output, "Candidates__both_lists_before_margin_filter.txt"))
  
  scs=screen_summary[,.(gene_id, origin, screen_DB, screen_SG, min_stat_overall,margin)]
  screen_1_selection=scs[margin>MARGIN]
  fwrite(screen_1_selection, paste0(folder_output, "Candidates__both_lists.txt"))
  return(screen_1_selection)
}


find_barcodes=function(screen_list, screen_dt,minBase=100,
                       COMP_MAIN_db, COMP_MAIN_sg, topN=10,
                       folder_output){
  
  # top N candidates
  screen_list=screen_list[1:topN]
  fwrite(screen_list, paste0(folder_output, "Selected_candidates_both_lists.txt"))
  
  # best from DB
  hits_db=screen_list[screen_DB=="DB"]$gene_id
  screen_dt_db=screen_dt[comparison==COMP_MAIN_db & gene_id %in% hits_db & baseMean>100]
  setorder(screen_dt_db, stat)
  screen_dt_db[, lapply(.SD, head, n=1), by=c("gene_id", "gene_name", "mutation_profile")]
  
  # best from SG
  hits_sg=screen_list[screen_SG=="SG"]$gene_id
  screen_dt_sg=screen_dt[comparison==COMP_MAIN_sg & gene_id %in% hits_sg & baseMean>100]
  setorder(screen_dt_sg, stat)
  screen_dt_sg=screen_dt_sg[, lapply(.SD, head, n=1), by=c("gene_id", "gene_name", "mutation_profile")]
  
  screen_barcodes=rbind(screen_dt_sg,screen_dt_db)
  setorder(screen_barcodes, stat)
  screen_barcodes[, lapply(.SD, head, n=1), by=c("gene_id", "gene_name", "mutation_profile")]
  fwrite(screen_barcodes, paste0(folder_output, "Best_barcodes_from_selected_candidates_both_lists.txt"))
  return(screen_barcodes)
}

plot_hits=function(screen_dt,COMP_MAIN_db, COMP_MAIN_sg, screen_list, folder_output){
  
  screen_dt=screen_dt[baseMean>100,]
  screen_dt_db=screen_dt[comparison==COMP_MAIN_db]
  screen_dt_sg=screen_dt[comparison==COMP_MAIN_sg]
  
  screen_paired=merge(screen_dt_db, screen_dt_sg,
                      by=c("barcode", "gene_id", "mutation_profile"),
                      suffixes=c("_db", "_sg"))
  
  for( GENE_SHOW in unique(screen_list$gene_id)){
    p_comp=ggplot(screen_paired[baseMean_db>100 | baseMean_sg>100], aes(stat_sg, stat_db))+
      geom_point(size = .5, alpha=.7, col="grey80")+
      geom_point(data=screen_paired[gene_id == GENE_SHOW], alpha=.7, size=1.5,aes(col=mutation_profile))+
      scale_color_manual(values=c(alt="red3", ref="black", "TAA"="gold"))+
      labs(y="Stat: Cluster TCRs vs control", x="Stat: Singlet TCRs vs control")+
      theme(legend.position = "none", axis.title = element_text(size=14))+
      geom_hline(yintercept = 0, lty="dashed", col="grey")+
      geom_vline(xintercept = 0, lty="dashed", col="grey")+
      ggtitle(paste0(GENE_SHOW, " in clusters vs singlets"))
    
    p_tot_db=ggplot(screen_dt_db, aes(baseMean, log2FoldChange))+
      geom_point( size = .5, alpha=.7, col="grey80")+
      scale_x_log10()+
      theme(legend.position = "none", axis.title = element_text(size=14))+
      geom_point(data=screen_dt_db[gene_id == GENE_SHOW], alpha=.7,size=1.5, aes(col=mutation_profile))+
      scale_color_manual(values=c(alt="red3", ref="black", "TAA"="gold"))+
      geom_hline(yintercept = 0, lty="dashed", col="grey")+
      ggtitle(paste0(GENE_SHOW, " in ", COMP_MAIN_db))+
      geom_text_repel(data=screen_dt_db[log2FoldChange<(-1)], aes(label = gene_id), size=2)
    
    p_tot_sg=ggplot(screen_dt_sg, aes(baseMean, log2FoldChange))+
      geom_point( size = .5, alpha=.7, col="grey80")+
      scale_x_log10()+
      theme(legend.position = "none", axis.title = element_text(size=14))+
      geom_point(data=screen_dt_sg[gene_id == GENE_SHOW], alpha=.7, size=1.5,aes(col=mutation_profile))+
      scale_color_manual(values=c(alt="red3", ref="black", "TAA"="gold"))+
      geom_hline(yintercept = 0, lty="dashed", col="grey")+
      ggtitle(paste0(GENE_SHOW, " in ", COMP_MAIN_sg))+
      geom_text_repel(data=screen_dt_sg[log2FoldChange<(-1)], aes(label = gene_id), size=2)
    
    cowplot::plot_grid(plotlist = list(p_comp, p_tot_db, p_tot_sg), align = "hv", nrow=1)
    ggsave(paste0(folder_output, "viz_results_",GENE_SHOW,".pdf"), width=16, height = 5)
  }
  
  p_tot_db_all=ggplot(screen_dt_db, aes(baseMean, log2FoldChange))+
    geom_point( size = .5, alpha=.7, col="grey80")+
    scale_x_log10()+
    theme(legend.position = "none", axis.title = element_text(size=14))+
    geom_point(data=screen_dt_db[gene_id %in% unique(screen_list$gene_id)], alpha=.7,size=1.5, aes(col=mutation_profile))+
    scale_color_manual(values=c(alt="red3", ref="black", "TAA"="gold"))+
    geom_hline(yintercept = 0, lty="dashed", col="grey")+
    ggtitle(paste0("All hits over ", COMP_MAIN_db))+
    geom_text_repel(data=screen_dt_db[gene_id %in% unique(screen_list$gene_id)],
                    aes(label = gene_id), size=2, max.overlaps = 10)+
    geom_line(data=screen_dt_db[gene_id %in% unique(screen_list$gene_id)],
              aes(group = gene_id), lwd=.3, col="black", lty="dashed")
  
  p_tot_sg_all=ggplot(screen_dt_sg, aes(baseMean, log2FoldChange))+
    geom_point( size = .5, alpha=.7, col="grey80")+
    scale_x_log10()+
    theme(legend.position = "none", axis.title = element_text(size=14))+
    geom_point(data=screen_dt_sg[gene_id %in% unique(screen_list$gene_id)],
               alpha=.7, size=1.5,aes(col=mutation_profile))+
    scale_color_manual(values=c(alt="red3", ref="black", "TAA"="gold"))+
    geom_hline(yintercept = 0, lty="dashed", col="grey")+
    ggtitle(paste0("All hits over ", COMP_MAIN_sg))+
    geom_text_repel(data=screen_dt_sg[gene_id %in% unique(screen_list$gene_id)],
                    aes(label = gene_id), size=2, max.overlaps = 10)+
    geom_line(data=screen_dt_sg[gene_id %in% unique(screen_list$gene_id)],
                    aes(group = gene_id), lwd=.3, col="grey30", lty="dashed")
  
  
  cowplot::plot_grid(plotlist = list(p_tot_db_all, p_tot_sg_all), align = "hv", nrow=1)
  ggsave(paste0(folder_output, "viz_results_overview.pdf"), width=11, height = 5)
}

# run pipeline
run_pipeline=function(screen_dt,
                      COMP_MAIN_db,
                      COMP_MAIN_sg,
                      pat_folder,
                      topN=topN,
                      MARGIN=2,
                      minBase=100){
  
  folder_output_dir=paste0(output_dir, pat_folder,"/")
  if(!dir.exists(folder_output_dir)){dir.create(folder_output_dir)}
  
  screen_clusters=hit_prioritisation(screen_dt = screen_dt, COMP_MAIN = COMP_MAIN_db,
                                     folder_output=folder_output_dir,
                                     minBase = minBase)
  screen_singlets=hit_prioritisation(screen_dt = screen_dt, COMP_MAIN = COMP_MAIN_sg,
                                     folder_output=folder_output_dir,
                                     minBase = minBase)
  screen_list=create_list(screen_clusters, screen_singlets,
                          folder_output=folder_output_dir,
                          MARGIN=MARGIN)
  
  screen_list=screen_list[1:topN,]
  screen_barcodes=find_barcodes(screen_list = screen_list,
                                  screen_dt = screen_dt,
                                  minBase = minBase, 
                                  COMP_MAIN_db = COMP_MAIN_db, 
                                  COMP_MAIN_sg = COMP_MAIN_sg,
                                  topN = topN,
                                folder_output=folder_output_dir)
  plot_hits(screen_dt=screen_dt, 
            COMP_MAIN_db = COMP_MAIN_db,
            COMP_MAIN_sg = COMP_MAIN_sg,
            screen_list=screen_barcodes,
            folder_output=folder_output_dir)
}

run_pipeline(screen_dt = screen_1,
             COMP_MAIN_db = "Cluster vs X1D3", 
             COMP_MAIN_sg = "Single vs X1D3",
             pat_folder = "P1",
             topN = 15,
             minBase = 100)

run_pipeline(screen_dt = screen_1,
             COMP_MAIN_db = "Cluster vs UT", 
             COMP_MAIN_sg = "Single vs UT",
             pat_folder = "P1_vs_UT",
             topN = 15,
             minBase = 100)

run_pipeline(screen_dt = screen_2a,
             COMP_MAIN_db = "Cluster vs UT", 
             COMP_MAIN_sg = "Single vs UT",
             pat_folder = "P5",
             topN = 15,
             minBase = 100)