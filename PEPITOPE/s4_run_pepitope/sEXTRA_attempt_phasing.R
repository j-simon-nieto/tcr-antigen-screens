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

htc_files=list.files("./PEPITOPE/s1_run_sarek/results_reseq/variant_calling/haplotypecaller/",recursive = T, pattern="haplotypecaller.filtered.vcf.gz$", full.names = F)
sample_names_htc=str_split_i(htc_files,"\\/",1)
names(htc_files)=sample_names_htc

# # rna fusion
fusion_vcf_files=list.files("./PEPITOPE/s3_run_rnafusion/results_reseq/vcf") # to be updated
sample_names_fusion=str_split_i(fusion_vcf_files,"_",1)
names(fusion_vcf_files)=sample_names_fusion

# # rna expr (RERAN INTO ONE!)
counts = readr::read_tsv("PEPITOPE/s2_run_rnaseq/results_initial_reseq/star_salmon/salmon.merged.gene_counts.tsv")
tpm = readr::read_tsv("PEPITOPE/s2_run_rnaseq/results_initial_reseq/star_salmon/salmon.merged.gene_tpm.tsv")

# # annotations
ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(ens106) = "UCSC"

# # metadata
meta=fread("PEPITOPE/s4_run_pepitope/match_rna_dna_meta.csv")

# set filters
min_rna_counts=1

# set output dir
out_path="PEPITOPE/s4_run_pepitope/reports_new_rna/"
if(!dir.exists(out_path)){dir.create(out_path)}

# set pat sample and run pipeline 
#pat_sample=meta$ID[meta$sample!="germline"][1]
estimate_impact=function(pat_sample){
  DNA_id=meta[ID==pat_sample]$DNA
  RNA_id=meta[ID==pat_sample]$RNA
  cat(pat_sample,"; dna: ",DNA_id,"; rna:",RNA_id,"\n")
    
  # get RNA expression for sample
  counts= counts %>% dplyr::select(gene_id, gene_name, count=RNA_id)
  tpm= tpm %>% dplyr::select(gene_id, gene_name, count=RNA_id)
  rna_sample = merge(counts, tpm, by=c("gene_id","gene_name"))
    
  colnames(rna_sample)=c("gene_id","gene_name","count","tpm")
  message("All good extracting the RNA")
    # ============
    
  # get & process variants
  vr1 = readVcfAsVRanges(
    paste0("./PEPITOPE/s1_run_sarek/results_reseq/variant_calling/mutect2/", mutect2_files[DNA_id]),
    use.names = TRUE)
  
  table(vr1$GT, vr1$PGT)
  mcols(vr1)$phase_summary= ifelse( vr1$GT == "0|1", "phased_haplo_2", 
                                    ifelse(vr1$GT == "1|0", "phased_haplo_1",
                                           "unphased"))
    
  sample_name = as.character(unique(sampleNames(vr1), value = TRUE))
  sel_sample_name=sample_name[grep(DNA_id,sample_name)]
  vr1 = filter_variants(vr1,sample = sel_sample_name,min_cov = 2, min_af = 0.05,pass = TRUE)
  
  # --- get the germline equivalent
  # no confidently phased variants with mutect2/haplotypecaller
  germline_sample_name=stringr::str_split_i(sample_name[!grepl(DNA_id,sample_name)],"_",2)
  
  vr_germline = readVcfAsVRanges(
    paste0("./PEPITOPE/s1_run_sarek/results_reseq/variant_calling/haplotypecaller/", htc_files[germline_sample_name]),
    use.names = TRUE)
  
  mcols(vr_germline)$phase_summary= ifelse(vr_germline$GT == "0|1", "phased_haplo_2", 
                                    ifelse(vr_germline$GT == "1|0", "phased_haplo_1",
                                           "unphased"))
  
  vr_germline = filter_variants(vr_germline,sample = sample_name[!grepl(DNA_id,sample_name)],min_cov = 2, min_af = 0.05,pass = TRUE)
  perc_phased_gl= 100-100*(sum(vr_germline$phase_summary=="unphased")/ length(vr_germline$phase_summary))
  cat(paste("Perc. phased variants in germline object",perc_phased_gl, "%\n", sep =" "))
  # ----------------------
    
  # ---- extra wrangling to fix discordant contigs
  seqlevelsStyle(vr1) <- "UCSC"
  common_seqlevels <- intersect(seqlevels(vr1), seqlevels(ens106))
  vr1_filtered <- keepSeqlevels(vr1, common_seqlevels, pruning.mode="coarse")
  
  # # for germline
  seqlevelsStyle(vr_germline) <- "UCSC"
  common_seqlevels <- intersect(seqlevels(vr_germline), seqlevels(ens106))
  vr_germline_filtered <- keepSeqlevels(vr_germline, common_seqlevels, pruning.mode="coarse")
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
  
  perc_phased= 100-100*(sum(subs$phase_summary=="unphased")/ length(subs$phase_summary))
  cat(paste("Perc. phased variants in final object",perc_phased, "%\n", sep =" "))
  
  # 1/ reporting germline background mutations
  # find peptides with germline mutations
  gr_gl <- granges(vr_germline_filtered[nchar(vr_germline_filtered@ref, keepNA = T) == 1 &  # ignore indels -- increased complexity
                                          nchar(vr_germline_filtered@alt, keepNA = F) == 1],
                   use.mcols = T)
  
  subs_93nt <- resize(subs,width = 93, fix = "center")
  
  germline_nearby <- subsetByOverlaps(subs_93nt,gr_gl, type="any")
  # annotate peptides with a mutation in the background
  mcols(subs)$germline_mut_bg=ifelse(subs$mut_id %in% germline_nearby$mut_id, "GL_BG_MUT", "")
  percentage_of_affected_peptides= 100*(sum(subs$germline_mut_bg=="GL_BG_MUT")/ length(subs$mut_id))
  cat(paste("Perc. mutated ids affected by germline backgroud mutations",percentage_of_affected_peptides, "%","(",sum(subs$germline_mut_bg=="GL_BG_MUT"),")\n", sep =" "))
  
  # 2/ how to integrate the mutations into the final peptides
  # 1. Standardize the Windowing
  subs_93nt <- resize(subs, width = 93, fix = "center")
  genome(subs_93nt) <- "hg38"
  
  # 2. Extract Reference Sequence (clean background)
  ref_seq_93 <- getSeq(asm, subs_93nt)
  
  # 3. Identify ONLY the somatic sites that have nearby germline mutations
  overlaps <- findOverlaps(subs_93nt, gr_gl, type = "any")
  query_hits <- unique(queryHits(overlaps))
  
  # 4. Initialize Background Sequences
  bg_seqs <- ref_seq_93
  
  # 5. Only loop through the rows that actually have overlaps
  for (i in query_hits) {
    gl_idx <- subjectHits(overlaps)[queryHits(overlaps) == i]
    gl_vars <- gr_gl[gl_idx]
    
    # Extract the DNAString object from the DNAStringSet
    tmp_str <- bg_seqs[[i]]
    
    for (j in seq_along(gl_vars)) {
      rel_pos <- start(gl_vars[j]) - start(subs_93nt[i]) + 1
      
      # Check window boundaries and protect the somatic center (pos 47)
      if (rel_pos >= 1 && rel_pos <= 93 && rel_pos != 47) {
        # FIX: Convert the ALT to a DNAString explicitly
        subseq(tmp_str, rel_pos, rel_pos) <- DNAString(as.character(gl_vars$alt[j]))
      }
    }
    # Place the modified DNAString back into the DNAStringSet
    bg_seqs[i] <- DNAStringSet(tmp_str)
  }
  
  # 6. Assign back to your object
  subs_gl <- subs
  subs_gl$ref_nuc <- bg_seqs
  
  # 7. Create Alt Nucleotides (Background + Somatic Mutation)
  alt_seqs <- bg_seqs
  # The center position (47) is currently the REF base. 
  # We replace it with the Somatic ALT base.
  subseq(alt_seqs, 47, 47) <- DNAStringSet(subs$alt)
  
  subs_gl$alt_nuc <- alt_seqs
  somatic_alts <- subs$alt # This is your "T"
  
  # Inject the somatic mutation at the center (position 47)
  # We use DNAStringSet to ensure compatibility
  subseq(alt_seqs, 47, 47) <- DNAStringSet(somatic_alts)
  subs_gl$alt_nuc <- alt_seqs
  
  # TO - DOs
  # --> actually THE RESULTS ARE WRONG
  # -> translation is simply wrong: take CDS into account
  # -> find a way to include in make_peptide and annotate with the
  # ---> fact that it's germline adaptaed
  # ---> what is the germline mutation
  

}

log_file <- file.path(out_path, "pepitope_germline_background_pipeline.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
lapply(meta$ID[meta$sample!="germline"],estimate_impact)
