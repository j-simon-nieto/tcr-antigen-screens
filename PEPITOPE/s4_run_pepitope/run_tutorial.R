library(pepitope)

library(pepitope)

# SNPs and small indels, e.g. from 'sarek' nf-core pipeline
variant_vcf_file = system.file("my_variants.vcf", package="pepitope")

# Combined fusion VCF file, e.g. from 'rnafusion' nf-core pipeline
fusion_vcf_file = system.file("my_fusions.vcf", package="pepitope")
v