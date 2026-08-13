library(data.table)
library(readxl)
setwd("/DATA/j.simon/JV/ANTIGEN_calling")

# load resources
orthoprimer=read_xlsx("PEPITOPE/s5_create_libraries/resources_zm/orthoprimer_patient_list.xlsx")
barcode_dt=read_xlsx("PEPITOPE/s5_create_libraries/resources_zm/minigene_barcode_list_all.xlsx")
Bbsl_SITE_1="GAAGAC"
CUT_SITE_1="TGCAAGC"
STOP="TGA"
CUT_SITE_2="CACTATG"	
Bbsl_SITE_2="GTCTTC"

# update: remove TAA antigens
taa_big=fread("PEPITOPE/s5_create_libraries/resources_zm/Large_TAA2.0_Oligo library - shared antigens V3 (93nt minigenes)_CSV_FINAL_LIBRARY.csv", header = T, skip = 4)
setDT(barcode_dt)
barcode_dt[, is_taa:= ifelse(Barcode %in% taa_big$`Oligo barcode`, "dup", "sg")]
barcode_dt=barcode_dt[!barcode_dt$Barcode %in% taa_big$`Oligo barcode`,]

peptide_files_dir="/DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s4_run_pepitope/reports_new_rna/"

# create a library
create_library=function(pX_agg,pat,primer_row){
  pX_library=data.table(
    replicate='1',
    patient_id=pat,
    detected_in=pX_agg$sample,
    gene_name=pX_agg$gene_name,
    mut_id=pX_agg$mut_id,
    pep_id=pX_agg$pep_id,
    pep_type=pX_agg$pep_type,
    n_tiles=pX_agg$n_tiles,
    ORTHOPRIMER_FW=orthoprimer$`Orthoprimer forward`[primer_row],
    Bbsl_SITE_1="GAAGAC",
    CUT_SITE_1="TGCAAGC",
    MINIGENE=pX_agg$tiled,
    STOP="TGA",
    BARCODE="empty_for_now", 
    CUT_SITE_2="CACTATG",
    Bbsl_SITE_2="GTCTTC",
    ORTHOPRIMER_REV=orthoprimer$`Orthoprimer reverse`[primer_row],
    nt_peptide=pX_agg$nt,
    peptide=pX_agg$peptide
  )
  
  pX_library_alt=data.table(
    replicate='2',
    patient_id=pat,
    detected_in=pX_agg$sample,
    gene_name=pX_agg$gene_name,
    mut_id=pX_agg$mut_id,
    pep_id=pX_agg$pep_id,
    pep_type=pX_agg$pep_type,
    n_tiles=pX_agg$n_tiles,
    ORTHOPRIMER_FW=orthoprimer$`Orthoprimer forward`[primer_row],
    Bbsl_SITE_1="GAAGAC",
    CUT_SITE_1="TGCAAGC",
    MINIGENE=pX_agg$tiled,
    STOP="TGA",
    BARCODE= "empty_for_now",
    CUT_SITE_2="CACTATG",
    Bbsl_SITE_2="GTCTTC",
    ORTHOPRIMER_REV=orthoprimer$`Orthoprimer reverse`[primer_row],
    nt_peptide=pX_agg$nt,
    peptide=pX_agg$peptide
  )

  
  pX_library_agg=rbind(pX_library,pX_library_alt)
  return(pX_library_agg)
}


p1_agg=fread(paste0(peptide_files_dir,"AKB459_aggregated_peptides.csv"))
p2_agg=fread(paste0(peptide_files_dir,"SJ2_aggregated_peptides.csv"))
p3_agg=fread(paste0(peptide_files_dir,"SJ18_aggregated_peptides.csv"))
p4_agg=fread(paste0(peptide_files_dir,"SJ20_aggregated_peptides.csv"))
p5_agg=fread(paste0(peptide_files_dir,"SJ32_aggregated_peptides.csv"))


library_all_pat=rbindlist(list(
  p1_library=create_library(p1_agg, pat = "P1", primer_row = 1),
  p2_library=create_library(p2_agg, pat = "P2", primer_row = 2),
  p3_library=create_library(p3_agg, pat = "P3", primer_row = 3),
  p4_library=create_library(p4_agg, pat = "P4", primer_row = 4),
  p5_library=create_library(p5_agg, pat = "P5", primer_row = 5)
))

library_all_pat[, BARCODE:= barcode_dt[1:nrow(library_all_pat),]$Barcode]
library_all_pat[, CONCATENATED_SEQ:= paste0(ORTHOPRIMER_FW,Bbsl_SITE_1,CUT_SITE_1,
                                           MINIGENE,STOP,BARCODE,CUT_SITE_2,
                                           Bbsl_SITE_2,ORTHOPRIMER_REV )]
library_all_pat[, seq_length:= nchar(CONCATENATED_SEQ)]


table(duplicated(library_all_pat$BARCODE))
table(library_all_pat$patient_id)
barcode_mat=dcast(library_all_pat, BARCODE~patient_id+replicate, fill = 0)

pdf("PEPITOPE/s5_create_libraries/new_rna_library_barcodes_heatmap.pdf", height = 5, width = 3)
ComplexHeatmap::Heatmap(barcode_mat[,-1], cluster_columns = T, cluster_rows = T, 
                        col = circlize::colorRamp2(colors=c("white", "coral3"), breaks = c(0,200)))
dev.off()

library_all_pat[,variant_name:= make.unique(paste0(patient_id,"_",pep_id,"_",pep_type,"_",replicate))]
fwrite(library_all_pat, "PEPITOPE/s5_create_libraries/new_rna_library_all_patients_duplo.txt")

fwrite(library_all_pat[,.(variant_name, CONCATENATED_SEQ)], "PEPITOPE/s5_create_libraries/new_rna_library_all_patients_duplo_order.txt")

probs=library_all_pat[library_all_pat$seq_length<174,]
nrow(probs)
nrow(probs)/nrow(library_all_pat)*100 # <1% has a short length
table(duplicated(library_all_pat$variant_name))

table(library_all_pat$replicate, library_all_pat$patient_id)
table(library_all_pat$ORTHOPRIMER_FW, library_all_pat$patient_id)
table(library_all_pat$seq_length, library_all_pat$patient_id)
table(library_all_pat$replicate, library_all_pat$patient_id,library_all_pat$detected_in)
table(library_all_pat$STOP)
table(library_all_pat$Bbsl_SITE_1,library_all_pat$Bbsl_SITE_2)
table(library_all_pat$CUT_SITE_1,library_all_pat$CUT_SITE_2)
nrow(library_all_pat)-length(unique(library_all_pat$BARCODE))

# compare to old results
lib_old=fread("PEPITOPE/s5_create_libraries/library_all_patients_duplo.txt")
length(intersect(library_all_pat$pep_id, lib_old$pep_id))
length(unique(library_all_pat$pep_id))

length(intersect(library_all_pat$CONCATENATED_SEQ, lib_old$CONCATENATED_SEQ))
length(setdiff(library_all_pat$CONCATENATED_SEQ, lib_old$CONCATENATED_SEQ))
length(unique(library_all_pat$CONCATENATED_SEQ))

for(i in names(library_all_pat)){
  print(i)
  print(length(intersect(library_all_pat[[i]], lib_old[[i]])))
  print(length(unique(library_all_pat[[i]])))
}

