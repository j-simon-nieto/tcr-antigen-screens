library(pepitope)
library(stringr)
library(data.table)
library(dplyr)

# mutect2 variant calling
mutect2_files=list.files("./PEPITOPE/s1_run_sarek/results_reseq/variant_calling/mutect2",
                         recursive = T, pattern="mutect2.filtered.vcf.gz$", full.names = F)

sample_names_mutect=str_split_i(mutect2_files,"_",1)
mutect2_data=lapply(mutect2_files, function(x){fread(paste0("./PEPITOPE/s1_run_sarek/results_reseq/variant_calling/mutect2/", x))})
names(mutect2_data)=sample_names_mutect
names(mutect2_files)=sample_names_mutect

# RNA fusions
fusion_vcf_files=list.files("./PEPITOPE/s3_run_rnafusion/results/vcf")
sample_names_fusion=str_split_i(fusion_vcf_files,"_",1)
fusion_data=lapply(fusion_vcf_files, function(x){fread(paste0("./PEPITOPE/s3_run_rnafusion/results/vcf/", x))})
names(fusion_data)=sample_names_fusion
names(fusion_vcf_files)=sample_names_fusion

# Preparation
# # Selecting the right reference genome
ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
asm = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(ens106) = "UCSC"

# Adding RNA expression 
counts_1 = readr::read_tsv("PEPITOPE/s2_run_rnaseq/results_initial/star_salmon/salmon.merged.gene_counts.tsv") |>
  dplyr::select(gene_id, gene_name, count=names(fusion_data)[1])
tpm_1 = readr::read_tsv("PEPITOPE/s2_run_rnaseq/results_initial/star_salmon/salmon.merged.gene_tpm.tsv") |>
  dplyr::select(gene_id, gene_name, tpm=names(fusion_data)[1])

# counts_2 = readr::read_tsv("PEPITOPE/s2_run_rnaseq/results_initial_2/star_salmon/salmon.merged.gene_counts.tsv") |>
#   dplyr::select(gene_id, gene_name, count=SAMPLE)
# tpm_2 = readr::read_tsv("PEPITOPE/s2_run_rnaseq/results_initial_2/star_salmon/salmon.merged.gene_tpm.tsv.tsv") |>
#   dplyr::select(gene_id, gene_name, tpm=SAMPLE)

rna_sample = inner_join(counts_1, tpm_1)
colnames(rna_sample)[1:2]=c("GENE_ID", "gene_name")


# # # ___ above process have to be systematic for each sample: create RNA/DNA metadata

vr1 = readVcfAsVRanges(paste0("./PEPITOPE/s1_run_sarek/results_reseq/variant_calling/mutect2/",mutect2_files["CF50136"])) |>
  filter_variants( sample = "AKB459_CF50136",min_cov=2, min_af=0.05, pass=TRUE)

# ---- extra wrangling
seqlevelsStyle(vr1) <- "UCSC"
common_seqlevels <- intersect(seqlevels(vr1), seqlevels(ens106))
vr1_filtered <- keepSeqlevels(vr1, common_seqlevels, pruning.mode="coarse")
# ----------------------

# --- extra wrangling to include expression as a filter
# --- changing gene name nomenclature
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# # 
gene_mapping <- biomaRt::getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  mart = ensembl)

setDT(gene_mapping)
gene_mapping <- gene_mapping[hgnc_symbol != ""]
setDT(rna_sample)
expressed_genes_hgnc <- rna_sample[count >= 5]$GENE_ID
expressed_ensembl <- gene_mapping[hgnc_symbol %in% expressed_genes_hgnc]$ensembl_gene_id

ann = annotate_coding(vr1_filtered, ens106, asm)

subs = ann |>
  subset_context(15) |>
  subset(GENEID %in% expressed_ensembl)

# READ FUSION VARIANTS
vr2 = readVcfAsVRanges(paste0("./PEPITOPE/s3_run_rnafusion/results/vcf/",fusion_vcf_files[1])) |>
  filter_fusions(min_reads=2, min_split_reads=1, min_tools=1)

# ---- extra wrangling
seqlevelsStyle(vr2) <- "UCSC"
common_seqlevels <- intersect(seqlevels(vr2), seqlevels(ens106))
vr2_filtered <- keepSeqlevels(vr2, common_seqlevels, pruning.mode="coarse")
genome(vr2_filtered) <- "GRCh38"
genome(asm) <- "GRCh38"
# ----------------------

fus = annotate_fusions(vr2_filtered, ens106, asm) |>
  subset_context_fusion(15)

# create construct tables and make a report
tiled = make_peptides(subs, fus) |>
  pep_tile() |>
  remove_cutsite(BbsI="GAAGAC")

report = make_report(ann, subs, fus, tiled)
writexl::write_xlsx(report, "./reseq_TEST_my_variants.xlsx")
