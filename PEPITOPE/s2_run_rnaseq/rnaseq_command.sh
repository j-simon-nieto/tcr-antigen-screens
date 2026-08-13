cd /DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s2_run_rnaseq/

nextflow run nf-core/rnaseq -r dev \
  --input samplesheet.csv \
  --outdir results_initial_reseq \
  -profile singularity \
  --aligner star_salmon \
  --email j.simon@nki.nl \
  --multiqc-title rna_dec25 \
  --genome GRCh38 \
  -work-dir /DATA/j.simon/JV/ANTIGEN_calling/PEPITOPE/s2_run_rnaseq/