nextflow run nf-core/rnafusion \
  -profile singularity \
  -c /DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s3_run_rnafusion/reseq_config.config \
  --tools all --no_cosmic \
  --input /DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s3_run_rnafusion/samplesheet.csv \
  --genomes_base /DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s3_run_rnafusion/references \
  --outdir /DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s3_run_rnafusion/results_reseq/ 