# =============================================================================
# Hit selection: explicit gates, robustness tiers, 4 output files per patient
# =============================================================================
# GOAL: select which hits from each screen to move forward with.
#
# Criteria:
#   TAA         : significant outlier dropout, consistent across its barcodes
#   Neoantigen  : significant outlier dropout in alt, no dropout in the WT ref
#   Robustness  : replicate support (n_sig) x effect size (best_lfc)
#
# Every decision is a named boolean column on the unit table, so a call can be
# audited without re-reading the filter chain.
#
#   barcode level : g_padj, g_lfc, g_out  ->  g3_passed
#   unit level    : n_sig, best_lfc, consistent, veto_*
#   tier          : A / B / C / vetoed / not_called
#
# Unit = gene_id (BLOCK level for fusions). For neoantigens this is
# mutation-level ONLY if gene_id encodes the mutation. If one gene carries >1
# mutation under the same gene_id, pass unit="pep_id" or the grain is wrong and
# it will fail silently.
# =============================================================================

library(data.table)
library(ggplot2)
library(ggrepel)
theme_set(theme_classic())
set.seed(150799)

output_dir <- "PEPITOPE/s7_NEW_screen_results/Hit_selection/stratified_output/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

.min_safe <- function(x) if (length(x) && any(!is.na(x))) min(x, na.rm = TRUE) else NA_real_
.max_safe <- function(x) if (length(x) && any(!is.na(x))) max(x, na.rm = TRUE) else NA_real_

.check_keys <- function(dt) {
  need <- c("gene_id", "pep_id", "ref_key", "is_fusion", "mutation_profile")
  miss <- setdiff(need, names(dt))
  if (length(miss))
    stop(sprintf("results table is missing %s -- re-run the screen script", paste(miss, collapse = ", ")))
  dt
}

# A fusion alt has two WT partners (5' and 3') when the upstream script emits
# ref_key_3p / ref_key_5p; otherwise there is a single ref_key. This resolves
# whichever columns are present so nothing downstream has to know which.
.rk_of <- function(x, units, unit = "gene_id") {
  kc <- intersect(c("ref_key_3p", "ref_key_5p", "ref_key"), names(x))
  if (!length(kc)) stop(".rk_of: no ref_key column found")
  rk <- unique(unlist(x[get(unit) %in% units & mutation_profile != "ref", ..kc]))
  rk[!is.na(rk)]
}

# ---------------------------------------------------------------------------
# 1. barcode-level gates
# ---------------------------------------------------------------------------
add_gates <- function(dt, padj_max, lfc_pass, n_sd) {
  dt <- copy(dt)
  # outlier threshold: n_sd below the mean stat of that comparison, computed on
  # the abundance-filtered set so it is not dragged down by low-count noise
  dt[, thr := mean(stat, na.rm = TRUE) - n_sd * sd(stat, na.rm = TRUE), by = comparison]
  dt[, g_padj := !is.na(padj) & padj < padj_max]
  dt[, g_lfc  := !is.na(log2FoldChange) & log2FoldChange < lfc_pass]
  dt[, g_out  := !is.na(stat) & stat <= thr]
  dt[, g3_passed := g_padj & g_lfc & g_out]
  dt[]
}

# ---------------------------------------------------------------------------
# 2. unit-level rollup
# ---------------------------------------------------------------------------
# test rows = barcodes the call is made on (TAA, or alt for neoantigens)
# ref  rows = WT partners, summarised separately and joined on ref_key so that
#             every block of a fusion finds the shared reference
# NOTE: ref-only gene_ids (e.g. PLXNA2.fs) no longer appear as units of their
#       own -- they exist only as the ref side of the join.
roll_units <- function(dt, by) {
  
  refs <- dt[mutation_profile == "ref",
             .(n_bc_ref = .N, n_sig_ref = sum(g3_passed),
               ref_best_stat = .min_safe(stat), ref_best_lfc = .min_safe(log2FoldChange)),
             by = .(ref_key)]
  
  tests <- dt[mutation_profile != "ref",
              .(class      = if (all(mutation_profile == "TAA")) "TAA" else "neo",
                is_fusion  = any(is_fusion),
                n_bc_test  = .N,
                n_sig      = sum(g3_passed),
                best_lfc   = .min_safe(log2FoldChange),
                worst_lfc  = .max_safe(log2FoldChange),
                best_stat  = .min_safe(stat),
                worst_stat = .max_safe(stat),
                best_padj  = .min_safe(padj)),
              by = by]
  
  # one row per (unit, WT partner): 1 for SNVs, 2 for fusions, 0 for TAA
  map <- unique(rbindlist(list(
    dt[mutation_profile != "ref" & !is.na(ref_key_3p), c(by, "ref_key_3p"), with = FALSE],
    dt[mutation_profile != "ref" & !is.na(ref_key_5p), c(by, "ref_key_5p"), with = FALSE]
  ), use.names = FALSE))
  setnames(map, 2L, "ref_key")
  
  agg <- refs[map, on = "ref_key", nomatch = NULL][
    , .(n_bc_ref = sum(n_bc_ref), n_sig_ref = sum(n_sig_ref),
        ref_best_stat = .min_safe(ref_best_stat),
        ref_best_lfc  = .min_safe(ref_best_lfc),
        n_ref_partners = .N), by = by]
  
  out <- agg[tests, on = by]
  out[, has_ref   := !is.na(n_bc_ref) & n_bc_ref > 0]
  out[, n_bc_ref  := fcoalesce(n_bc_ref,  0L)]
  out[, n_sig_ref := fcoalesce(n_sig_ref, 0L)]
  out[, n_ref_partners := fcoalesce(n_ref_partners, 0L)]
  out[]
}

# ---------------------------------------------------------------------------
# 3. per-comparison prioritisation  (returns, does not write)
# ---------------------------------------------------------------------------
hit_prioritisation <- function(screen_dt, COMP_MAIN, COMP_CTRL = "UT vs B",
                               minBase = 30,
                               padj_max = 0.05,
                               lfc_pass = -0.2,
                               lfc_high = -0.5,
                               n_sd = 2,
                               unit = "gene_id") {
  
  stopifnot(COMP_MAIN %in% screen_dt$comparison, COMP_CTRL %in% screen_dt$comparison)
  
  dt    <- add_gates(.check_keys(screen_dt[baseMean > minBase]), padj_max, lfc_pass, n_sd)
  
  m <- roll_units(dt[comparison == COMP_MAIN], by = unit)
  k <- roll_units(dt[comparison == COMP_CTRL], by = unit)[
    , .(unit_ = get(unit), ctrl_n_sig = n_sig, ctrl_best_stat = best_stat)]
  setnames(k, "unit_", unit)
  m <- k[m, on = unit]
  
  # -- vetoes ---------------------------------------------------------------
  # neo: one ref barcode clearing all three gates.
  m[, veto_wt_ref := class == "neo" & n_sig_ref > 0]
  m[, veto_donor  := class == "TAA" & !is.na(ctrl_n_sig) & ctrl_n_sig > 0] # not enforced
  m[, veto_TAA_consistency := class == "TAA" & !is.na(worst_lfc) & worst_lfc > 0]
  # flagged, NOT vetoed: mutation-specificity is untestable without a WT partner,
  # so the ref test was skipped rather than passed
  m[, no_ref := class == "neo" & !has_ref]
  m[, veto        := veto_wt_ref | veto_TAA_consistency]
  m[, veto_reason := fcase(veto_wt_ref, "wt_ref_dropout",
                           veto_TAA_consistency, "inconsistent_dropout",
                           default = NA_character_)]
  
  # -- robustness axes ------------------------------------------------------
  m[, support_tier := fcase(n_sig >= 2, "replicated",
                            n_sig == 1, "single",
                            default    = "none")]
  m[, effect_tier  := fcase(is.na(best_lfc),      NA_character_,
                            best_lfc <  lfc_high, "high",
                            best_lfc <  lfc_pass, "moderate",
                            default              = "weak")]
  m[, consistent := !is.na(worst_stat) & worst_stat < 0]
  
  m[, tier := fcase(
    veto,                                                     "vetoed",
    n_sig == 0,                                               "not_called",
    support_tier == "replicated" & effect_tier == "high",     "A",
    support_tier == "replicated" & effect_tier == "moderate", "B",
    support_tier == "single"     & effect_tier == "high",     "B",
    support_tier == "single"     & effect_tier == "moderate", "C",
    default = "not_called")]
  m[, tier := factor(tier, levels = c("A","B","C","vetoed","not_called"), ordered = TRUE)]
  
  # margin: neo -> vs its WT partner (via ref_key); TAA -> vs the donor control.
  # NA (not Inf) when the partner is absent -- see no_ref.
  m[, margin        := fifelse(class == "neo", ref_best_stat  - best_stat,
                               ctrl_best_stat - best_stat)]
  m[, margin_source := fifelse(class == "neo", "ref", "donor_ctrl")]
  
  setorder(m, tier, best_stat)
  message(sprintf("[%s]\n%s", COMP_MAIN,
                  paste(capture.output(print(m[, .N, by = .(class, tier)][order(class, tier)])), collapse = "\n")))
  m[]
}

# ---------------------------------------------------------------------------
# 4. combine the two TCR pools
# ---------------------------------------------------------------------------
create_list <- function(screen_clusters, screen_singlets, unit = "gene_id") {
  s <- merge(screen_clusters, screen_singlets,
             by = c(unit, "class"), suffixes = c("_DB", "_SG"), all = TRUE)
  
  lv <- levels(s$tier_DB)
  s[, best_tier := factor(lv[pmin(as.integer(tier_DB), as.integer(tier_SG), na.rm = TRUE)],
                          levels = lv, ordered = TRUE)]
  s[, n_pools_called := (!is.na(tier_DB) & tier_DB %in% c("A","B","C")) +
      (!is.na(tier_SG) & tier_SG %in% c("A","B","C"))]
  s[, best_stat_overall := pmin(best_stat_DB, best_stat_SG, na.rm = TRUE)]
  s[, best_lfc_overall  := pmin(best_lfc_DB,  best_lfc_SG,  na.rm = TRUE)]
  s[, margin_overall    := pmax(margin_DB,    margin_SG,    na.rm = TRUE)]
  # fusion parent, so blocks of one fusion can be capped or grouped downstream
  s[, fusion_parent := fifelse(grepl("--", get(unit), fixed = TRUE),
                               sub("-[0-9]+$", "", get(unit)), NA_character_)]
  
  setorder(s, best_tier, -n_pools_called, best_stat_overall)
  s[]
}

# ---------------------------------------------------------------------------
# 5. one representative barcode per selected unit
# ---------------------------------------------------------------------------
find_barcodes <- function(selected, gates_dt, COMP_MAIN_db, COMP_MAIN_sg,
                          unit = "gene_id") {
  
  rk <- .rk_of(gates_dt, selected[[unit]], unit)
  
  bc <- gates_dt[comparison %in% c(COMP_MAIN_db, COMP_MAIN_sg) &
                   (get(unit) %in% selected[[unit]] |
                      (mutation_profile == "ref" & ref_key %in% rk))]
  setorder(bc, stat)
  test <- bc[mutation_profile != "ref"][, .SD[1], by = unit]      # one per unit
  ref  <- bc[mutation_profile == "ref"][, .SD[1], by = "ref_key"] # one per WT partner
  bc   <- rbind(test, ref, use.names = TRUE, fill = TRUE)
  setorder(bc, stat)
  bc[]
}

# ---------------------------------------------------------------------------
# 6. plots
# ---------------------------------------------------------------------------
plot_hits <- function(gates_dt, gp_dt, selected, COMP_MAIN_db, COMP_MAIN_sg,
                      COMP_CTRL = "UT vs B", minBase = 30,
                      unit = "gene_id", folder_output) {
  
  pal <- c(alt = "red3", ref = "black", "TAA" = "gold")
  d   <- gates_dt[baseMean > minBase]
  d_db <- d[comparison == COMP_MAIN_db]
  d_sg <- d[comparison == COMP_MAIN_sg]
  d_ct <- d[comparison == COMP_CTRL]
  
  join_cols <- c("barcode", unit, "mutation_profile",
                 intersect(c("ref_key", "ref_key_3p", "ref_key_5p", "is_fusion"),
                           names(d_db)))
  paired <- merge(d_db, d_sg, by = join_cols, suffixes = c("_db", "_sg"))
  
  # highlight the unit's own barcodes plus every WT barcode controlling it
  hi_rows <- function(x, g, rk) x[get(unit) == g |
                                    (mutation_profile == "ref" & ref_key %in% rk)]
  
  ma_panel <- function(x, ttl, sel_unit, sel_rk) {
    hi <- hi_rows(x, sel_unit, sel_rk)
    ggplot(x, aes(baseMean, log2FoldChange)) +
      geom_point(size = .5, alpha = .7, col = "grey80") + scale_x_log10() +
      geom_point(data = hi, aes(col = mutation_profile), size = 1.5, alpha = .7) +
      scale_color_manual(values = pal) +
      geom_hline(yintercept = 0, lty = "dashed", col = "grey") +
      geom_hline(yintercept = -0.2, lwd = .3, lty = "dashed", col = "grey") +
      geom_text_repel(data = hi, aes(label = .data[[unit]]), size = 2.5) +
      theme(legend.position = "none", axis.title = element_text(size = 12)) +
      ggtitle(paste0(sel_unit, " | ", ttl))
  }
  
  per_gene <- file.path(folder_output, "per_unit_plots")
  tier_of <- setNames(as.character(selected$best_tier), selected[[unit]])
  rank_of <- setNames(sprintf("%02d", seq_len(nrow(selected))), selected[[unit]])
  for (tt in unique(tier_of))
    dir.create(file.path(per_gene, paste0("tier_", tt)), recursive = TRUE, showWarnings = FALSE)
  
  for (g in unique(selected[[unit]])) {
    
    rk <- .rk_of(d, g, unit)          # both WT partners for a fusion, one otherwise
    
    p1 <- ggplot(paired, aes(stat_sg, stat_db)) +
      geom_point(size = .5, alpha = .7, col = "grey80") +
      geom_point(data = hi_rows(paired, g, rk), aes(col = mutation_profile),
                 size = 1.5, alpha = .7) +
      scale_color_manual(values = pal) +
      geom_hline(yintercept = 0, lty = "dashed", col = "grey") +
      geom_vline(xintercept = 0, lty = "dashed", col = "grey") +
      theme(legend.position = "none", axis.title = element_text(size = 12)) +
      labs(y = paste0("Stat: ", COMP_MAIN_db), x = paste0("Stat: ", COMP_MAIN_sg)) +
      ggtitle(paste0(g, " | pools"))
    
    # gp_dt carries ref_key; matching a WT by gene_id stopped working once
    # gene_id became variant-level
    p5 <- ggplot(hi_rows(gp_dt, g, rk), aes(variable, norm + 1)) +
      geom_point(aes(col = mutation_profile), size = 4, alpha = .5) +
      geom_line(aes(group = guide)) + scale_y_log10() +
      scale_color_manual(values = c(alt = "red3", ref = "grey", "TAA" = "gold")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
            legend.position = "none") +
      labs(x = NULL, y = "normalised counts + 1") + ggtitle(paste0(g, " | profile"))
    
    cowplot::plot_grid(p1,
                       ma_panel(d_ct, COMP_CTRL,    g, rk),
                       ma_panel(d_db, COMP_MAIN_db, g, rk),
                       ma_panel(d_sg, COMP_MAIN_sg, g, rk),
                       p5, nrow = 1)
    ggsave(file.path(per_gene, paste0("tier_", tier_of[[g]]),
                     paste0("viz_", gsub("[^A-Za-z0-9._-]", "_", g),
                            "__rank_", rank_of[[g]], ".pdf")),
           width = 20, height = 4)
  }
  
  overview <- function(x, ttl) {
    sel <- x[get(unit) %in% selected[[unit]]]
    ggplot(x, aes(baseMean, log2FoldChange)) +
      geom_point(size = .5, alpha = .3, col = "grey80") + scale_x_log10() +
      geom_point(data = sel, aes(col = mutation_profile), size = 1.5, alpha = .7) +
      geom_line(data = sel, aes(group = .data[[unit]]), lwd = .3,
                lty = "dashed", col = "grey30") +
      geom_text_repel(data = sel, aes(label = .data[[unit]]), size = 2, max.overlaps = 10) +
      scale_color_manual(values = pal) +
      geom_hline(yintercept = c(0, -0.2), lty = "dashed", col = "grey", lwd = c(.5, .3)) +
      theme(legend.position = "none", axis.title = element_text(size = 12)) +
      ggtitle(paste0("All selected | ", ttl))
  }
  cowplot::plot_grid(overview(d_db, COMP_MAIN_db), overview(d_sg, COMP_MAIN_sg),
                     align = "hv", nrow = 1)
  ggsave(file.path(folder_output, "viz_overview.pdf"), width = 11, height = 5)
}

# ---------------------------------------------------------------------------
# 7. pipeline  -- writes exactly 4 tables
# ---------------------------------------------------------------------------
run_pipeline <- function(screen_dt, COMP_MAIN_db, COMP_MAIN_sg, pat_folder,
                         gp_dt, topN = 30, minBase = 30, minBase_plot = 30,
                         COMP_CTRL = "UT vs B", padj_max = 0.05,
                         lfc_pass = -0.2, lfc_high = -0.5, n_sd = 2,
                         unit = "gene_id") {
  
  fo <- file.path(output_dir, pat_folder); dir.create(fo, recursive = TRUE, showWarnings = FALSE)
  
  gates <- add_gates(.check_keys(screen_dt[baseMean > minBase]), padj_max, lfc_pass, n_sd)
  
  args <- list(screen_dt = screen_dt, COMP_CTRL = COMP_CTRL, minBase = minBase,
               padj_max = padj_max, lfc_pass = lfc_pass, lfc_high = lfc_high,
               n_sd = n_sd, unit = unit)
  db <- do.call(hit_prioritisation, c(args, COMP_MAIN = COMP_MAIN_db))
  sg <- do.call(hit_prioritisation, c(args, COMP_MAIN = COMP_MAIN_sg))
  
  cand    <- create_list(db, sg, unit = unit)
  called  <- cand[best_tier %in% c("A","B","C")]
  selected <- head(called, topN)
  bcs     <- find_barcodes(selected, gates, COMP_MAIN_db, COMP_MAIN_sg, unit = unit)
  
  # ---- 4 output tables ----------------------------------------------------
  # 1: full audit trail, every unit, both pools
  fwrite(cand, file.path(fo, "01_units_all.txt"))
  # 2: the decision table -- tiered hits only
  fwrite(called, file.path(fo, "02_candidates_tiered.txt"))
  # 3: top N with the barcode to carry forward
  fwrite(merge(selected, bcs[, .(get(unit), barcode, gene_name, mutation_profile,
                                 comparison, baseMean, log2FoldChange, stat, padj)],
               by.x = unit, by.y = "V1", all.x = TRUE, sort = FALSE),
         file.path(fo, "03_selected_with_barcode.txt"))
  # 4: barcode-level evidence behind the selection, including the WT refs that
  #    the selected units are controlled by (matched via ref_key)
  sel_rk <- .rk_of(gates, selected[[unit]], unit)
  fwrite(gates[(get(unit) %in% selected[[unit]] |
                  (mutation_profile == "ref" & ref_key %in% sel_rk)) &
                  comparison %in% c(COMP_MAIN_db, COMP_MAIN_sg, COMP_CTRL)],
                  file.path(fo, "04_barcode_evidence.txt"))
  
  plot_hits(gates, gp_dt, selected, COMP_MAIN_db, COMP_MAIN_sg,
            COMP_CTRL = COMP_CTRL, minBase = minBase_plot, unit = unit,
            folder_output = fo)
  
  message(sprintf("[%s] %d units -> %d called -> %d selected -> %d barcodes | no_ref in called: %d",
                  pat_folder, nrow(cand), nrow(called), nrow(selected), nrow(bcs),
                  called[no_ref_DB %in% TRUE | no_ref_SG %in% TRUE, .N]))
  invisible(list(units = cand, called = called, selected = selected, barcodes = bcs))
}

# ---------------------------------------------------------------------------
# 8. run
# ---------------------------------------------------------------------------
screen_1  <- fread("PEPITOPE/s7_NEW_screen_results/output_screen_1/updated_screen_results/screen1_results.txt")
screen_2a <- fread("PEPITOPE/s7_NEW_screen_results/output_screen_2a/updated_screen_results/screen2a_results.txt")

run_pipeline(screen_dt    = screen_1,
             COMP_MAIN_db = "Cluster vs UT",
             COMP_MAIN_sg = "Single vs UT",
             pat_folder   = "P1",
             topN         = 30,
             minBase      = 100,
             gp_dt = fread("PEPITOPE/s7_NEW_screen_results/output_screen_1/updated_screen_results/gene_profiles.txt"))

run_pipeline(screen_dt    = screen_2a,
             COMP_MAIN_db = "Cluster vs UT",
             COMP_MAIN_sg = "Single vs UT",
             pat_folder   = "P5",
             topN         = 30,
             minBase      = 100,
             gp_dt = fread("PEPITOPE/s7_NEW_screen_results/output_screen_2a/updated_screen_results/gene_profiles.txt"))

run_pipeline(screen_dt    = screen_1,
             COMP_MAIN_db = "X1D3 vs UT",
             COMP_MAIN_sg = "X1D3 vs UT",
             pat_folder   = "P1_1d3_control",
             topN         = 30,
             minBase      = 100,
             gp_dt = fread("PEPITOPE/s7_NEW_screen_results/output_screen_1/updated_screen_results/gene_profiles.txt"))