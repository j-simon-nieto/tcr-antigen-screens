# PENDING QUESTIONS
#
# MultiQC highlights duplication quality issues, issue? -> no
# Variant calling is based only on Mutect2?
# Is it ok to update the annotation references?
# Must rnafusion be run with cosmic?

library(pepitope)
library(stringr)
library(data.table)
library(dplyr)
library(ggplot2)

#================ VARIANT CALLING ================

# load requirements
# # variant calling
mutect2_files=list.files("./PEPITOPE/s1_run_sarek/results_reseq/variant_calling/mutect2",recursive = T, pattern="mutect2.filtered.vcf.gz$", full.names = F)
sample_names_mutect=str_split_i(mutect2_files,"_",1)
names(mutect2_files)=sample_names_mutect

# # rna fusion
fusion_vcf_files=list.files("./PEPITOPE/s3_run_rnafusion/results/vcf")
sample_names_fusion=str_split_i(fusion_vcf_files,"_",1)
names(fusion_vcf_files)=sample_names_fusion

# # rna expr (RERUN INTO ONE!)
counts_1 = readr::read_tsv("PEPITOPE/s2_run_rnaseq/results_initial/star_salmon/salmon.merged.gene_counts.tsv")
tpm_1 = readr::read_tsv("PEPITOPE/s2_run_rnaseq/results_initial/star_salmon/salmon.merged.gene_tpm.tsv")
counts_2 = readr::read_tsv("PEPITOPE/s2_run_rnaseq/results_initial_2/star_salmon/salmon.merged.gene_counts.tsv")
tpm_2 = readr::read_tsv("PEPITOPE/s2_run_rnaseq/results_initial_2/star_salmon/salmon.merged.gene_tpm.tsv")

# # annotations
ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(ens106) = "UCSC"

# # metadata
meta=fread("PEPITOPE/s4_run_pepitope/match_rna_dna_meta.csv")

# set filters
min_rna_counts=1

# set output dir
out_path="PEPITOPE/s4_run_pepitope/reports_mincount1_reseq/"
if(!dir.exists(out_path)){dir.create(out_path)}

# create function
make_pep_report=function(pat_sample){
  DNA_id=meta[ID==pat_sample]$DNA
  RNA_id=meta[ID==pat_sample]$RNA
  message(pat_sample,"; dna: ",DNA_id,"; rna:",RNA_id)
  
  # get RNA expression for sample
  if (RNA_id %in% colnames(counts_1)){
    counts_1=counts_1 %>% dplyr::select(gene_id, gene_name, count=RNA_id)
    tpm_1=tpm_1 %>% dplyr::select(gene_id, gene_name, count=RNA_id)
    rna_sample = merge(counts_1, tpm_1, by=c("gene_id","gene_name"))
  } else {
    counts_2=counts_2 %>% dplyr::select(gene_id, gene_name, count=RNA_id)
    tpm_2=tpm_2 %>% dplyr::select(gene_id, gene_name, count=RNA_id)
    rna_sample = merge(counts_2, tpm_2, by=c("gene_id","gene_name"))
  }
  
  colnames(rna_sample)=c("gene_id","gene_name","count","tpm")
  message("All good extracting the RNA")
  # ============
  
  # get & process variants
  vr1 = readVcfAsVRanges(
    paste0("./PEPITOPE/s1_run_sarek/results_reseq/variant_calling/mutect2/", mutect2_files[DNA_id]),
    use.names = TRUE)
  
  sample_name = as.character(unique(sampleNames(vr1), value = TRUE))
  sel_sample_name=sample_name[grep(DNA_id,sample_name)]
  vr1 = filter_variants(vr1,sample = sel_sample_name,min_cov = 2, min_af = 0.05,pass = TRUE)
  
  # ---- extra wrangling to fix discordant contigs
  seqlevelsStyle(vr1) <- "UCSC"
  common_seqlevels <- intersect(seqlevels(vr1), seqlevels(ens106))
  vr1_filtered <- keepSeqlevels(vr1, common_seqlevels, pruning.mode="coarse")
  # ----------------------
  
  # ---- filter by expression
  setDT(rna_sample)
  expressed_genes_hgnc <- rna_sample[count >= min_rna_counts]$gene_id
  # ----------------------
  
  ann = annotate_coding(vr1_filtered, ens106, asm)
  subs = ann |>
    subset(gene_name %in% expressed_genes_hgnc) |>
    subset_context(15) 
  message("All good extracting the variants")
  # ============
  
  # get & process FUSIONS
  vr2 = readVcfAsVRanges(paste0("./PEPITOPE/s3_run_rnafusion/results/vcf/",fusion_vcf_files[RNA_id])) |>
    filter_fusions(min_reads=min_rna_counts, min_split_reads=1, min_tools=1)
  
  # ---- extra wrangling
  #seqlevelsStyle(vr2) <- "UCSC"
  common_seqlevels <- intersect(seqlevels(vr2), seqlevels(ens106))
  vr2_filtered <- keepSeqlevels(vr2, common_seqlevels, pruning.mode="coarse")
  genome(vr2_filtered) <- "GRCh38"
  genome(asm) <- "GRCh38"
  # ----------------------
  
  fus = annotate_fusions(vr2_filtered, ens106, asm) |>
    subset_context_fusion(15)
  message("All good extracting the fusions")
  # ============
  
  # tile cDNA and create the report
  tiled = make_peptides(subs, fus) |>
    pep_tile() |>
    remove_cutsite(BbsI="GAAGAC")
  
  report = make_report(ann, subs, fus, tiled)
  fwrite(report$`93 nt Peptides`, paste0(out_path,pat_sample,"_peptides.txt"))
  writexl::write_xlsx(report, paste0(out_path,pat_sample,"_report_file.xlsx"))
  message(pat_sample, "--REPORTED\n")
}

lapply(meta$ID[meta$sample!="germline"],make_pep_report)


### mini recap
res=list.files("/DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s4_run_pepitope/reports_mincount1_reseq/",
           full.names = T, pattern = "peptides.txt")
names_res=gsub("_peptides.txt","",str_split_i(res,"//",2))
res_list=lapply(res, function(x){message(x);fread(x)})
names(res_list)=names_res

res_dt=rbindlist(res_list, idcol = "sample")
res_dt[, pat:= str_split_i(sample, "_",1)]
res_dt[, origin:= gsub("_NA","",paste(str_split_i(sample, "_",2),
                        str_split_i(sample, "_",3),
                        str_split_i(sample, "_",4),sep="_"))]
numbers=res_dt[,.N,by=c("sample","pep_type","pat","origin")]

ggplot(numbers[pep_type=="alt"], aes(origin, N, fill=origin))+
  geom_col(position="dodge", width=.8)+theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
  scale_fill_manual(values=c("cell_line"="lightyellow2","cell_line_P21"="lightyellow2",
                             "cell_line_P6"="lightyellow2","cell_line_P19"="lightyellow2",
                             "cell_line_P4"='lightyellow2',"tumor"="lightblue3"))+
  facet_wrap(~pat, scales="free_x", nrow=1, space="free_x")
ggsave(paste0(out_path,"tot_alt_peptides_report_file.pdf"), width=5, height = 3)

mut_ids_list=lapply(res_list, function(x){unique(x$mut_id)})

pdf(paste0(out_path,"intersect_peptides_report_file.pdf"), width=8, height = 10)
cowplot::plot_grid(
  ggVennDiagram::ggVennDiagram(mut_ids_list[grep("AK", names(mut_ids_list))], label_alpha = 0)+scale_fill_gradient(low = "white", high = "coral3")+
    theme(legend.position = "none")+ scale_x_continuous(expand = expansion(mult = .4)),
  ggVennDiagram::ggVennDiagram(mut_ids_list[grep("SJ2_", names(mut_ids_list))], label_alpha = 0)+scale_fill_gradient(low = "white", high = "coral3")+
    theme(legend.position = "none",
          text = element_text(size = 4))+ scale_x_continuous(expand = expansion(mult = .4)),
  ggVennDiagram::ggVennDiagram(mut_ids_list[grep("SJ20", names(mut_ids_list))], label_alpha = 0)+scale_fill_gradient(low = "white", high = "coral3")+
    theme(legend.position = "none",
          text = element_text(size = 4))+ scale_x_continuous(expand = expansion(mult = .4)),
  ggVennDiagram::ggVennDiagram(mut_ids_list[grep("SJ32", names(mut_ids_list))], label_alpha = 0)+scale_fill_gradient(low = "white", high = "coral3")+
    theme(legend.position = "none",
          text = element_text(size = 4))+ scale_x_continuous(expand = expansion(mult = .4)),
  ggVennDiagram::ggVennDiagram(mut_ids_list[grep("SJ18", names(mut_ids_list))], label_alpha = 0)+scale_fill_gradient(low = "white", high = "coral3")+
    theme(legend.position = "none",
          text = element_text(size = 4))+ scale_x_continuous(expand = expansion(mult = .4)),
  nrow=3
)
dev.off()



#================ AGGREGATE PER PATIENT ================

# AKB459 ______
AKB459_cell_line_specific=setdiff(res_list$AKB459_cell_line$pep_id, res_list$AKB459_tumor$pep_id)
AKB459_tumor_specific=setdiff( res_list$AKB459_tumor$pep_id,res_list$AKB459_cell_line$pep_id)

AKB459_agg=rbindlist(idcol = "sample", list("Cell_line"=res_list$AKB459_cell_line[pep_id %in% AKB459_cell_line_specific],
                                 "Tumor"=res_list$AKB459_tumor))

AKB459_agg[, sample:= ifelse(sample=="Tumor" & !(pep_id %in% AKB459_tumor_specific), "shared", sample)]
# # # checks (pep id and mut id slightly different because of number of peptides?)
table(AKB459_agg$sample)
AKB459_agg[, length(unique(mut_id)), by="sample"]

fwrite(AKB459_agg, "AKB459_aggregated_peptides.csv")

# SJ2 ______
SJ2_cell_lineP4_vs_tumor=setdiff(res_list$SJ2_cell_line_P4$pep_id, res_list$SJ2_tumor$pep_id)

SJ2_P4_tum_agg=rbindlist(idcol = F, list("Cell_line"=res_list$SJ2_cell_line_P4[pep_id %in% SJ2_cell_lineP4_vs_tumor],
                                            "Tumor"=res_list$SJ2_tumor))

SJ2_cell_lineP19_vs_others=setdiff(res_list$SJ2_cell_line_P19$pep_id,SJ2_P4_tum_agg$pep_id)

SJ2_agg=rbindlist(idcol = "sample", list("Cell_line_P19"=res_list$SJ2_cell_line_P19[pep_id %in% SJ2_cell_lineP19_vs_others],
                                                "shared"=SJ2_P4_tum_agg))
SJ2_cell_lineP4_vs_others=setdiff(res_list$SJ2_cell_line_P4$pep_id, union(res_list$SJ2_tumor$pep_id,res_list$SJ2_cell_line_P19$pep_id))
SJ2_tumor_vs_others=setdiff(res_list$SJ2_tumor$pep_id, union(res_list$SJ2_cell_line_P4$pep_id,res_list$SJ2_cell_line_P19$pep_id))

SJ2_agg[, sample:= ifelse(sample=="shared" & pep_id %in% SJ2_cell_lineP4_vs_others, "Cell_line_P4",
                          ifelse(sample=="shared" & pep_id %in% SJ2_tumor_vs_others, "Tumor",sample))]
# # # checks (pep id and mut id slightly different because of number of peptides?)
table(SJ2_agg$sample)
SJ2_agg[, length(unique(mut_id)), by="sample"]

fwrite(SJ2_agg, "SJ2_aggregated_peptides.csv")

# SJ18 ______
SJ18_cell_lineP6_vs_tumor=setdiff(res_list$SJ18_cell_line_P6$pep_id, res_list$SJ18_tumor$pep_id)

SJ18_P6_tum_agg=rbindlist(idcol = F, list("Cell_line"=res_list$SJ18_cell_line_P6[pep_id %in% SJ18_cell_lineP6_vs_tumor],
                                         "Tumor"=res_list$SJ18_tumor))

SJ18_cell_lineP21_vs_others=setdiff(res_list$SJ18_cell_line_P21$pep_id,SJ18_P6_tum_agg$pep_id)

SJ18_agg=rbindlist(idcol = "sample", list("Cell_line_P21"=res_list$SJ18_cell_line_P21[pep_id %in% SJ18_cell_lineP21_vs_others],
                                         "shared"=SJ18_P6_tum_agg))
SJ18_cell_lineP6_vs_others=setdiff(res_list$SJ18_cell_line_P6$pep_id, union(res_list$SJ18_tumor$pep_id,res_list$SJ18_cell_line_P21$pep_id))
SJ18_tumor_vs_others=setdiff(res_list$SJ18_tumor$pep_id, union(res_list$SJ18_cell_line_P6$pep_id,res_list$SJ18_cell_line_P21$pep_id))

SJ18_agg[, sample:= ifelse(sample=="shared" & pep_id %in% SJ18_cell_lineP6_vs_others, "Cell_line_P6",
                          ifelse(sample=="shared" & pep_id %in% SJ18_tumor_vs_others, "Tumor",sample))]
# # # checks (pep id and mut id slightly different because of number of peptides?)
table(SJ18_agg$sample)
SJ18_agg[, length(unique(mut_id)), by="sample"]

fwrite(SJ18_agg, "SJ18_aggregated_peptides.csv")

# SJ20 ______
SJ20_cell_line_specific=setdiff(res_list$SJ20_cell_line$pep_id, res_list$SJ20_tumor$pep_id)
SJ20_tumor_specific=setdiff( res_list$SJ20_tumor$pep_id,res_list$SJ20_cell_line$pep_id)

SJ20_agg=rbindlist(idcol = "sample", list("Cell_line"=res_list$SJ20_cell_line[pep_id %in% SJ20_cell_line_specific],
                                            "Tumor"=res_list$SJ20_tumor))

SJ20_agg[, sample:= ifelse(sample=="Tumor" & !(pep_id %in% SJ20_tumor_specific), "shared", sample)]
# # # checks (pep id and mut id slightly different because of number of peptides?)
table(SJ20_agg$sample)
SJ20_agg[, length(unique(mut_id)), by="sample"]

fwrite(SJ20_agg, "SJ20_aggregated_peptides.csv")

# SJ32 ______
SJ32_cell_line_specific=setdiff(res_list$SJ32_cell_line$pep_id, res_list$SJ32_tumor$pep_id)
SJ32_tumor_specific=setdiff( res_list$SJ32_tumor$pep_id,res_list$SJ32_cell_line$pep_id)

SJ32_agg=rbindlist(idcol = "sample", list("Cell_line"=res_list$SJ32_cell_line[pep_id %in% SJ32_cell_line_specific],
                                          "Tumor"=res_list$SJ32_tumor))

SJ32_agg[, sample:= ifelse(sample=="Tumor" & !(pep_id %in% SJ32_tumor_specific), "shared", sample)]
# # # checks (pep id and mut id slightly different because of number of peptides?)
table(SJ32_agg$sample)
SJ32_agg[, length(unique(mut_id)), by="sample"]

fwrite(SJ32_agg, "SJ32_aggregated_peptides.csv")

#================ QUALITY CONTROL ================