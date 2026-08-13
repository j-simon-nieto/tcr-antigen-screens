# NOTES:
# the TAA libraries have this one barcode that takes half of the reads
# this is not a used barcode
# the barcode is created as a combination of the remaining cut site and the annexed primer
# apparently the reads are assigned to that barcode instead of the intended one

# HYPOTHESIS
# two barcodes are being called from one single read
# this is because the amount of good reads called are the same as the templates detected
# also, why removing the barcode doesnt impact the results
# also the depth matches the expected one if taking into account that the amount of minigenes is x6
# --> the one thing we can't explain is why this second barcode is detected despite being in the wrong position
# --> also doubting why we only find it in the TAAs --> something to do with the barcodes included?

# EXTRA
# word file explains how to interpret the fasta reads
# --> start with the new annexxed primer and end with the 7 first nucleotides of the minigene

# TO DO
# run without R to enforce exact matching
# check if the "bad" barcode is detected 3M only if allowing the mismatch
# --> this is because it depends on a G--> C in the remaining cut site

# CURRENT SOLUTION:
# before enforcing the exact matching I trimmed the inital part of the sequences until the 26th
# --> barcode + stop + minigene
# this works but we will try to bypass it with the manual matching

# FINAL RESOLUTION:
# manual run worked, our hypothesis was correct


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

demux_metrics=fread(paste0(temp_dir, "/demux-metrics.txt"))
dmm=melt(demux_metrics, id.vars=c("sample_id", "barcode"))
ggplot(dmm, aes(sample_id, value))+
  facet_wrap(~variable, scale="free_y", nrow=1)+
  geom_col(fill="coral3")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
ggsave("PEPITOPE/s6_library_runs/first_run_march_26/output/sanity_metrics.pdf",
       height=4, width = 10)

# 2/ Counting construct barcodes
# ================================================

# # Need all valid barcodes
library(readxl)
barcode_dt=read_xlsx("PEPITOPE/s5_create_libraries/resources_zm/minigene_barcode_list_all.xlsx")
early_barcodes=barcode_dt$Barcode

# # Need all constructs is a sample divided list

# # # 2.2.1 patient libraries
all_constructs=fread("PEPITOPE/s5_create_libraries/new_rna_library_all_patients_duplo.txt")
setnames(all_constructs, "MINIGENE", "tiled")
setnames(all_constructs, "BARCODE", "barcode")
const_list=split(all_constructs, f=all_constructs$patient_id)

# # # TAA lists
taa_big=fread("PEPITOPE/s6_library_runs/shared_by_JV/large_TAA.csv")
taa_small=fread("PEPITOPE/s6_library_runs/shared_by_JV/small_TAA_fix.csv")

# *rename big list
setnames(taa_big, "BARCODE", "barcode")
setnames(taa_big, "GENE NAME 2", "mut_id")
setnames(taa_big, "MINIGENE SEQ", "tiled")
taa_big[, gene_name:= stringr::str_split_i(mut_id, "-", 1)]
taa_big[, pep_id:= mut_id]
taa_big[, pep_type:= "TAA"]
taa_big[, sample_id:= "TAA_big"]
taa_big[, patient_id:= "TAA_big"]
taa_big=taa_big[!is.na(LENGTH)]

# *rename small list
setnames(taa_small, "Oligo barcode", "barcode")
setnames(taa_small, "GENE NAME", "gene_name")
setnames(taa_small, "MINIGENE SEQ", "tiled")
taa_small[, mut_id:=paste0(gene_name, "_", `OLIGO NR`)]
taa_small[, pep_id:= mut_id]
taa_small[, pep_type:= "TAA"]
taa_small[, sample_id:= "TAA_small"]
taa_small[, patient_id:= "TAA_small"]
# exact repetition of the same row? 4764 EBV_EBNA1   ACTGTCA
taa_small=taa_small[!duplicated(taa_small)]

const_list[["TAA_big"]]=taa_big
const_list[["TAA_small"]]=taa_small

# we have duplicate barcodes so we need to make custom code
# this is because there's overlap across patients and TAA libraries
# will not be a problem for the readouts
# dset = count_bc(tdir = temp_dir, 
#                 all_constructs=const_list,
#                 valid_barcodes=valid_barcodes)

custom_count_external=function(tdir, valid_barcodes, reverse_complement) {
  tsv = tibble(name = as.character(valid_barcodes)) |>
    mutate(barcode=name, gene=name)
  if (reverse_complement)
    tsv$barcode = as.character(reverseComplement(DNAStringSet(tsv$barcode)))
  lpath = file.path(tdir, "lib.tsv")
  utils::write.table(tsv, file=lpath, sep="\t", row.names=FALSE, quote=FALSE)
  
  fqs = list.files(tdir, pattern="\\.R1\\.fq\\.gz$", full.names=TRUE)
  if (length(fqs) == 0)
    stop("No fastq files found to count in directory ", sQuote(tdir))
  
  guideCounterWrapper::guidecounter_count(
    input = fqs,
    library = lpath,
    offset_min_fraction = 0.2,
    output = file.path(tdir, "barcodes")
  )
  
  res = c(counts="barcodes.counts.txt", stats="barcodes.stats.txt") |>
    lapply(\(f) file.path(tdir, f)) |>
    lapply(readr::read_tsv, show_col_types=FALSE)
}

custom_count_bc=function(tdir, all_constructs, valid_barcodes,
                         reverse_complement=TRUE, sample,
                         barcode_start) { # the first posisition is the 26th
  samples = readr::read_tsv(file.path(tdir, "samples.tsv"), show_col_types=FALSE)
  samples=samples[samples$sample_id %in% sample,]
  all_samples = strsplit(samples$patient, "+", fixed=TRUE) |> unlist()
  missing = setdiff(all_samples, names(all_constructs))
  if (length(missing) > 0)
    stop("Missing minigene annotations for: ", paste(sQuote(missing), collapse=", "))
  
  construct_df = pepitope:::merge_constructs(all_constructs)
  if (missing(valid_barcodes))
    valid_barcodes = construct_df$barcode
  if (!is.character(valid_barcodes) && !is.factor(valid_barcodes))
    stop("'valid_barcodes' must be a character vector")
  if (!all(construct_df$barcode %in% valid_barcodes))
    stop("'all_constructs' contains barcodes not in 'valid_barcodes'")
  
  # remove the beggining of the reads to avoid stumbling on the barcode that is created
  # Trim FASTQs to start at barcode position
  fqs <- list.files(tdir, pattern="\\.R1\\.fq\\.gz$", full.names=TRUE)
  trimmed_dir <- file.path(tdir, "trimmed")
  dir.create(trimmed_dir, showWarnings=FALSE)
  
  lapply(fqs, function(fq) {
    out <- file.path(trimmed_dir, basename(fq))
    cmd <- sprintf(
      "zcat %s | awk 'NR%%4==2||NR%%4==0 {$0=substr($0,%d)} {print}' | gzip > %s",
      fq, barcode_start, out)
    system(cmd)
  })
  
  # copy samples.tsv to trimmed_dir so count_external finds it
  file.copy(file.path(tdir, "samples.tsv"), trimmed_dir)
  file.copy(file.path(tdir, "lib.tsv"),     trimmed_dir)  # if pre-existing
  # ______
  
  res = custom_count_external(trimmed_dir, valid_barcodes, reverse_complement) # requires installation of guide-counter in conda env
  
  stats = res$stats |> mutate(sample_id = sub("\\.R1$", "", label))
  counts = data.matrix(res$counts[-(1:2)])
  rownames(counts) = res$counts$guide
  colnames(counts) = sub("\\.R1$", "", colnames(counts))
  
  meta = samples |>
    left_join(stats |> select(sample_id, total_reads, mapped_reads)) |>
    mutate(patient = factor(patient),
           smp = paste(ifelse(nchar(origin)>8, stringr::word(origin, 1), origin), rep, sep="-"),
           short = paste(patient, smp),
           label = sprintf("%s (%s)", short, sample_id))
  
  rows = data.frame(barcode=valid_barcodes) |>
    left_join(construct_df) |>
    mutate(bc_type = ifelse(is.na(bc_type), "unused", bc_type),
           bc_type = factor(bc_type, levels=c(names(all_constructs), "unused")))
  
  SummarizedExperiment(
    list(counts = counts[, meta$sample_id, drop=FALSE]),
    colData = meta,
    rowData = rows
  )
}

ext_valid_barcodes=union(early_barcodes, union(taa_big$barcode, taa_small$barcode))

res_count_bc=lapply(unique(sample_sheet$sample_id), function(x){
  const_list=const_list[names(const_list)==x]
  res=custom_count_bc(tdir = temp_dir, 
                      all_constructs=const_list,
                      valid_barcodes=ext_valid_barcodes,
                      sample=x)
  return(res)
})


if(!dir.exists("PEPITOPE/s6_library_runs/first_run_march_26/tmp_files/")){
  dir.create("PEPITOPE/s6_library_runs/first_run_march_26/tmp_files/")}
names(res_count_bc)=unique(sample_sheet$sample_id)
R.utils::copyDirectory( from=temp_dir, 
                        to="PEPITOPE/s6_library_runs/first_run_march_26/tmp_files/", 
                        recursive=T)

saveRDS(res_count_bc, "PEPITOPE/s6_library_runs/first_run_march_26/output/per_sample_counts.rds")

# 3/ Quality Control plots
# ==============================================

# # # Overview
library(plotly)
plot_list=list(p1=plot_reads(res_count_bc$P1),
     p2=plot_reads(res_count_bc$P2),
     p3=plot_reads(res_count_bc$P3),
     p4=plot_reads(res_count_bc$P4),
     p5=plot_reads(res_count_bc$P5),
     TAA_big=plot_reads(res_count_bc$TAA_big),
     TAA_small=plot_reads(res_count_bc$TAA_small))

if(!dir.exists("PEPITOPE/s6_library_runs/first_run_march_26/output/magic/")){
  dir.create("PEPITOPE/s6_library_runs/first_run_march_26/output/magic/")}

pdf("PEPITOPE/s6_library_runs/first_run_march_26/output/read_dsitribution_per_patient.pdf",
    height=12, width = 14)
cowplot::plot_grid(plotlist = plot_list, align = "hv", nrow = 4)
dev.off()

lapply(names(res_count_bc), function(x){
  dir_out="PEPITOPE/s6_library_runs/first_run_march_26/output/magic/"
  plot_dir= paste0(dir_out,x,"_distr_reads.html")
  y=ggplotly(plot_distr(res_count_bc[[x]]), height=500, tooltip="text")
  htmlwidgets::saveWidget(y, plot_dir, selfcontained = T)
})


list_of_distr=lapply(names(res_count_bc), function(x){
  y=plot_distr(res_count_bc[[x]])+theme_classic()
  return(y)
})

pdf("PEPITOPE/s6_library_runs/first_run_march_26/output/barcode_abundance_distribution_per_patient.pdf",
    height=12, width = 8)
cowplot::plot_grid(plotlist = list_of_distr, align = "hv", nrow = 4)
dev.off()

mat_list=lapply(res_count_bc, function(x){
  y=as.data.table(assay(x), keep.rownames = "barcode")
  y[,patient_id:= names(y)[2]]
  names(y)[2]="counts"
  
  z=as.data.table(rowData(x))
  w=merge(y,z, by=c("barcode"))
  return(w)
  })

mat_dt=rbindlist(mat_list, fill=T, idcol = "library_id")

mat_dt[, bc_type_2:= ifelse(bc_type=="unused", "unused", "selected")]
mat_dt[, tot_counts:=sum(counts), by="library_id"]
mat_dt[, rel_counts:=counts/tot_counts]
fwrite(mat_dt, "PEPITOPE/s6_library_runs/first_run_march_26/output/metadata_n_counts.txt")

# # barcode counts distribution
ggplot(mat_dt, aes(library_id, counts, fill = bc_type_2))+
  geom_boxplot()+scale_y_log10()+theme_classic()+
  ggrepel::geom_text_repel(data=mat_dt[counts>3e6],
            aes(label=barcode), size=3, min.segment.length = 0.01)
ggsave("PEPITOPE/s6_library_runs/first_run_march_26/output/counts_boxplots_outlier.pdf",
    height=5, width = 7)

mat_dt[, median(counts), by=c("library_id", "bc_type_2")]

# # proportion of reads
mat_dt[, count_lvl:= cut(counts, breaks = c(0,1,10,100,1000,Inf), include.lowest = T, ordered_result = T)]
ggplot(mat_dt, aes(bc_type_2,fill=count_lvl))+geom_bar(position='fill')+
  theme_classic()+facet_wrap(~library_id, nrow=1)+
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=.5))+
  scale_fill_manual(values = c("red4", "red1", "orange", "yellow", "palegreen3"))
ggsave("PEPITOPE/s6_library_runs/first_run_march_26/output/counts_size_distribution.pdf",
       height=5, width = 7)


# # inspecting the size of the contamination
ggplot(mat_dt, aes(bc_type_2, y=rel_counts,fill = barcode=="TGACTACGTTGA"))+
  geom_col(position="stack")+
  theme_classic()+facet_wrap(~library_id, nrow=1)+
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=.5))
ggsave("PEPITOPE/s6_library_runs/first_run_march_26/output/total_reads_outlier_highlight.pdf",
       height=5, width = 7)

mat_dt[counts>100000]
mat_dt[barcode=="TGACTACGTTGA"]



# 4. Additional analysis to understand what's going on with the extra barcode
# ============================================

# # Trying to discard that it is a cross contamination/ issue with our code modifications
#====
# KEEP GROUPS:
barcodes_list=lapply(const_list, function(x){x$barcode})
ggVennDiagram::ggVennDiagram(barcodes_list, show_intersect = T) # indeed non overlapping

# # all patients
patient_ids=c("P1", "P2", "P3", "P4", "P5")
const_list_pats=const_list[names(const_list) %in% patient_ids]
dset_patients=custom_count_bc(tdir = temp_dir, 
                      all_constructs=const_list_pats,
                      valid_barcodes=ext_valid_barcodes,
                      sample=patient_ids)
plot_reads(dset_patients)
plot_distr(dset_patients)+
  geom_hline(yintercept = c(10,100,1000), lty="dashed", col="grey")+
  geom_vline(xintercept = c(0.05), lty="dashed", col="grey")+
  theme_classic()+scale_y_continuous(breaks=c(0,10,100,1000,10000), transform = "log10")+
  scale_x_continuous(breaks=c(0,0.05,0.5,1))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

# # patients + TAA_big
Pat_bigTAA_ids=c("TAA_big", patient_ids)
const_list_bigTAA=const_list[names(const_list) %in% Pat_bigTAA_ids]
dset_TAA_big=custom_count_bc(tdir = temp_dir, 
                              all_constructs=const_list_bigTAA,
                              valid_barcodes=ext_valid_barcodes,
                              sample=Pat_bigTAA_ids, barcode_start = 26)
plot_reads(dset_TAA_big)
plot_distr(dset_TAA_big)+
  geom_hline(yintercept = c(10,100,1000), lty="dashed", col="grey")+
  geom_vline(xintercept = c(0.05), lty="dashed", col="grey")+
  theme_classic()+scale_y_continuous(breaks=c(0,10,100,1000,10000), transform = "log10")+
  scale_x_continuous(breaks=c(0,0.05,0.5,1))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

ggplotly(plot_distr(dset_TAA_big), height=500, tooltip="text")

# # # manual checks here
counts=as.data.table(assay(dset_TAA_big), keep.rownames = "barcode")
meta_row=as.data.table(rowData(dset_TAA_big))
counts_row=merge(meta_row, counts, by="barcode")


# # patients + TAA_small
Pat_smallTAA_ids=c("TAA_small", patient_ids)
const_list_smallTAA=const_list[names(const_list) %in% Pat_smallTAA_ids]
dset_TAA_small=custom_count_bc(tdir = temp_dir, 
                             all_constructs=const_list_smallTAA,
                             valid_barcodes=ext_valid_barcodes,
                             sample=Pat_smallTAA_ids)
plot_reads(dset_TAA_small)
plot_distr(dset_TAA_small)+
  geom_hline(yintercept = c(10,100,1000), lty="dashed", col="grey")+
  geom_vline(xintercept = c(0.05), lty="dashed", col="grey")+
  theme_classic()+scale_y_continuous(breaks=c(0,10,100,1000,10000), transform = "log10")+
  scale_x_continuous(breaks=c(0,0.05,0.5,1))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

ggplotly(plot_distr(dset_TAA_small), height=500, tooltip="text")
# the barcode still is found only in TAA libraries and takes most reads.

# # What happens if we don't include the barcode?
exempt_dset_TAA_small=custom_count_bc(tdir = temp_dir, 
                               all_constructs=const_list_smallTAA,
                               valid_barcodes=setdiff(ext_valid_barcodes, "TGACTACGTTGA"),
                               sample=Pat_smallTAA_ids)
plot_reads(exempt_dset_TAA_small)
plot_distr(exempt_dset_TAA_small)+
  geom_hline(yintercept = c(10,100,1000), lty="dashed", col="grey")+
  geom_vline(xintercept = c(0.05), lty="dashed", col="grey")+
  theme_classic()+scale_y_continuous(breaks=c(0,10,100,1000,10000), transform = "log10")+
  scale_x_continuous(breaks=c(0,0.05,0.5,1))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
ggplotly(plot_distr(exempt_dset_TAA_small), height=500, tooltip="text")
# We observe the same thing: (no outlier) but much lower read than in the other libraries


# # EXPLORE IF IT'S A CONTAMINATION
library(Biostrings)
exp_barcode="TGACTACGTTGA"
other_seqs=list(Bbsl_SITE_1="GAAGAC",
CUT_SITE_1="TGCAAGC",
STOP="TGA",
CUT_SITE_2="CACTATG",
Bbsl_SITE_2="GTCTTC")
RC_barcode=Biostrings::reverseComplement(Biostrings::DNAString(exp_barcode))
random_barcode="AAGCGGCAACGT"
random_RC_barcode=Biostrings::reverseComplement(Biostrings::DNAString(random_barcode))

ctr_barcode2="GTAGCTGAAGTT"
RC_ctr_barcode2=Biostrings::reverseComplement(Biostrings::DNAString(ctr_barcode2))

check="AAGAGATTCGAG"
Biostrings::reverseComplement(Biostrings::DNAString(check))

RC_list=lapply(other_seqs, function(x){
  RC_barcode=Biostrings::reverseComplement(Biostrings::DNAString(x))
})

# taa big primer:
fw_primer_big=unique(taa_big$`FW PRIMER`)
RC_fw_primer_big=Biostrings::reverseComplement(Biostrings::DNAString(fw_primer_big))

rev_primer_big=unique(taa_big$`REV PRIMER`)
RC_rev_primer_big=Biostrings::reverseComplement(Biostrings::DNAString(rev_primer_big))


# BASH --> paste
# for sample in TAA_big TAA_small; do
# zcat /DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s6_library_runs/first_run_march_26/tmp_files/${sample}.R1.fq.gz \
# | awk 'NR%4==1 {header=$0} NR%4==2 {seq=$0}
#            NR%4==2 && (seq ~ /AACTTCAGCTAC/) {
#              print header "\n" seq
#            }' \
# > /DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s6_library_runs/first_run_march_26/output/${sample}_ctr2_reads.fa
# done
