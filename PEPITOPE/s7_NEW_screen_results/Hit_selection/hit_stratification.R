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
#
# -----------------------------------------------------------------------------
# 2026-10-09 fixes (a-f):
#   (a) roll_units now emits `ag_class` (neo | TAA); viral antigens fold into
#       TAA and are kept as candidates, judged under the same conditions.
#       All downstream `class=="neo"` / bare `mutation_profile` refs that
#       crashed hit_prioritisation now use ag_class.
#   (b) TAA WT key is transported into ref_key_3p at load (.transport_taa_refkey)
#       so the existing map/roll_units machinery pairs TAA alts with their WT.
#   (c) create_list merges on `unit` only, then coalesces class / ag_class, so a
#       unit whose class flips between pools is no longer split into half-NA rows.
#   (d) best_tier = pmin across pools kept intentionally (best hit across the
#       board; a veto in one pool is NOT terminal).
#   (e) .check_keys now also requires ref_key_3p / ref_key_5p (indexed
#       unconditionally by roll_units).
#   (f) non-fatal load-time check that TAA ref_keys survived upstream integration.
#   margin is now ref-based whenever a WT partner exists (neo OR TAA), donor
#   control only as fallback -- consistent with treating TAA/viral like neo.
# =============================================================================

library(data.table)
library(ggplot2)
library(ggrepel)
theme_set(theme_classic())
set.seed(150799)

output_dir <- "PEPITOPE/s7_NEW_screen_results/Hit_selection/stratified_output_TAAwt/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

.min_safe <- function(x) if (length(x) && any(!is.na(x))) min(x, na.rm = TRUE) else NA_real_
.max_safe <- function(x) if (length(x) && any(!is.na(x))) max(x, na.rm = TRUE) else NA_real_

# (e) ref_key_3p / ref_key_5p added: roll_units indexes them unconditionally,
#     so a results table that predates the fusion work should fail here loudly.
.check_keys <- function(dt) {
  need <- c("gene_id", "pep_id", "ref_key", "ref_key_3p", "ref_key_5p",
            "is_fusion", "mutation_profile")
  miss <- setdiff(need, names(dt))
  if (length(miss))
    stop(sprintf("results table is missing %s -- re-run the screen script", paste(miss, collapse = ", ")))
  dt
}

# (b) TAA/viral test rows carry their WT key in `ref_key` but NOT in
#     ref_key_3p/5p (those are NA for non-fusions and were never populated for
#     TAAs). Move it into ref_key_3p so roll_units' `map` -- which only reads
#     ref_key_3p/5p -- finds the WT partner. This is the minimal transport you
#     asked for; nothing downstream has to learn about TAAs.
.transport_taa_refkey <- function(dt) {
  dt <- copy(dt)
  if (!"ref_key_3p" %in% names(dt)) dt[, ref_key_3p := NA_character_]
  if (!"ref_key_5p" %in% names(dt)) dt[, ref_key_5p := NA_character_]
  dt[grepl("TAA|Viral", mutation_profile) & !grepl("ref", mutation_profile) &
       is.na(ref_key_3p) & is.na(ref_key_5p) & !is.na(ref_key),
     ref_key_3p := ref_key]
  dt[]
}

# (f) non-fatal guard: if the results predate the TAA ref_key/type_TAA
#     integration, every TAA has ref_key = NA and all TAA control silently
#     vanishes. Warn rather than stop (the run is still valid for neoantigens).
.check_taa_integration <- function(dt, tag) {
  n_taa <- dt[grepl("TAA|Viral", mutation_profile) & !grepl("ref", mutation_profile), .N]
  n_key <- dt[grepl("TAA|Viral", mutation_profile) & !grepl("ref", mutation_profile) &
                !is.na(ref_key), .N]
  if (n_taa > 0 && n_key == 0)
    warning(sprintf("[%s] %d TAA/viral test rows but 0 carry a ref_key -- results likely predate the TAA integration; TAA WT control will be absent.",
                    tag, n_taa))
  else
    message(sprintf("[%s] TAA ref_key check: %d/%d TAA/viral test rows keyed to a WT.",
                    tag, n_key, n_taa))
  invisible(dt)
}

source("PEPITOPE/s7_NEW_screen_results/palette_screens.R")

# A fusion alt has two WT partners (5' and 3') when the upstream script emits
# ref_key_3p / ref_key_5p; otherwise there is a single ref_key. This resolves
# whichever columns are present so nothing downstream has to know which.
.rk_of <- function(x, units, unit = "gene_id") {
  kc <- intersect(c("ref_key_3p", "ref_key_5p", "ref_key"), names(x))
  if (!length(kc)) stop(".rk_of: no ref_key column found")
  rk <- unique(unlist(x[get(unit) %in% units & !grepl("ref", mutation_profile), ..kc]))
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
# test rows = barcodes the call is made on (TAA/viral, or alt for neoantigens)
# ref  rows = WT partners, summarised separately and joined on ref_key so that
#             every block of a fusion finds the shared reference
# NOTE: ref-only gene_ids (e.g. PLXNA2.fs) no longer appear as units of their
#       own -- they exist only as the ref side of the join.
roll_units <- function(dt, by) {
  
  refs <- dt[grepl("ref", mutation_profile),
             .(n_bc_ref = .N, n_sig_ref = sum(g3_passed),
               ref_best_stat = .min_safe(stat), ref_best_lfc = .min_safe(log2FoldChange)),
             by = .(ref_key)]
  
  # (a) ag_class: neo vs TAA (viral folds into TAA). One row per unit, so all()
  #     over the non-ref test rows of that unit is well defined.
  tests <- dt[!grepl("ref", mutation_profile),
              .(ag_class   = if (all(grepl("TAA|Viral", mutation_profile))) "TAA" else "neo",
                class      = if (all(!is.na(ref_key))) "controlled" else "no_ref",
                is_fusion  = any(is_fusion),
                n_bc_test  = .N,
                n_sig      = sum(g3_passed),
                best_lfc   = .min_safe(log2FoldChange),
                worst_lfc  = .max_safe(log2FoldChange),
                best_stat  = .min_safe(stat),
                worst_stat = .max_safe(stat),
                best_padj  = .min_safe(padj)),
              by = by]
  
  # one row per (unit, WT partner): 1 for SNVs, 2 for fusions, 1 for keyed TAAs
  # (ref_key_3p now carries the TAA WT key via .transport_taa_refkey), 0 otherwise
  map <- unique(rbindlist(list(
    dt[!grepl("ref", mutation_profile) & !is.na(ref_key_3p), c(by, "ref_key_3p"), with = FALSE],
    dt[!grepl("ref", mutation_profile) & !is.na(ref_key_5p), c(by, "ref_key_5p"), with = FALSE]
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
  # (a) all antigen-type logic keyed on ag_class / has_ref, not the dead
  #     class=="neo" literal.
  # neo or keyed TAA: a WT barcode of this unit clears all three gates.
  m[, veto_wt_ref := has_ref & n_sig_ref > 0]
  # computed, NOT enforced (reported for assessment, per instruction)
  m[, veto_donor  := ag_class == "TAA" & !is.na(ctrl_n_sig) & ctrl_n_sig > 0]
  m[, veto_TAA_consistency := ag_class == "TAA" & !is.na(worst_lfc) & worst_lfc > 0]
  # flagged, NOT vetoed: mutation-specificity untestable without a WT partner
  m[, no_ref := ag_class == "neo" & !has_ref]
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
  
  # margin: vs the WT partner whenever one exists (neo OR keyed TAA); vs the
  # donor control only as fallback. NA (not Inf) when neither is available.
  m[, margin        := fifelse(has_ref, ref_best_stat  - best_stat,
                               ctrl_best_stat - best_stat)]
  m[, margin_source := fifelse(has_ref, "ref", "donor_ctrl")]
  
  setorder(m, tier, best_stat)
  message(sprintf("[%s]\n%s", COMP_MAIN,
                  paste(capture.output(print(m[, .N, by = .(ag_class, tier)][order(ag_class, tier)])), collapse = "\n")))
  m[]
}

# ---------------------------------------------------------------------------
# 4. combine the two TCR pools
# ---------------------------------------------------------------------------
# (c) merge on `unit` only; class / ag_class are unit invariants, so coalesce
#     them rather than merging on a column that can flip between pools and split
#     a unit into two half-NA rows.
create_list <- function(screen_clusters, screen_singlets, unit = "gene_id") {
  s <- merge(screen_clusters, screen_singlets,
             by = unit, suffixes = c("_DB", "_SG"), all = TRUE)
  s[, class    := fcoalesce(class_DB,    class_SG)]
  s[, ag_class := fcoalesce(ag_class_DB, ag_class_SG)]
  
  # (d) best across pools kept intentionally: a veto in one pool does not sink a
  #     unit that is a clean hit in the other.
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
  
  s[, rank_stat := frank(best_stat_overall, ties.method = "min")]
  s[, rank_lfc  := frank(best_lfc_overall,  ties.method = "min")]
  s[, rank_best := pmin(rank_stat, rank_lfc)]
  setorder(s, rank_best, best_stat_overall)
  s[, sel_pos := .I]   # now the post-sort rank (moved after setorder)
  s[]
}

# ---------------------------------------------------------------------------
# 5. one representative barcode per selected unit
# ---------------------------------------------------------------------------
find_barcodes <- function(selected, gates_dt, COMP_MAIN_db, COMP_MAIN_sg,
                          unit = "gene_id") {
  
  comps <- c(COMP_MAIN_db, COMP_MAIN_sg)
  
  # -- 1. strongest alt/test barcode per unit, pooled across DB + Single ------
  # ranked purely on log2FoldChange (most negative = strongest dropout). The
  # winning row carries the comparison it won in, which the ref must match.
  tst <- gates_dt[comparison %in% comps &
                    get(unit) %in% selected[[unit]] &
                    !grepl("ref", mutation_profile) &
                    !is.na(log2FoldChange)]
  setorder(tst, log2FoldChange)
  test <- tst[, .SD[1], by = unit]                 # winner + its comparison
  
  # annotation only (neo | TAA); the ref rule below is uniform
  if ("ag_class" %in% names(selected))
    test <- selected[, c(unit, "ag_class"), with = FALSE][test, on = unit]
  
  # -- 2. matched WT barcodes: for EACH WT partner of the unit, the best ref ---
  #    barcode from the SAME comparison the alt won in. Fusions -> 5' and 3';
  #    TAAs with a WT -> their reference; partnerless units -> none.
  ref_list <- lapply(seq_len(nrow(test)), function(i) {
    u    <- test[[unit]][i]
    cmp  <- test$comparison[i]
    rk_u <- .rk_of(gates_dt, u, unit)              # all WT key(s) for this unit
    if (!length(rk_u)) return(NULL)
    r <- gates_dt[comparison == cmp & grepl("ref", mutation_profile) &
                    ref_key %in% rk_u & !is.na(log2FoldChange)]
    if (!nrow(r)) return(NULL)
    setorder(r, log2FoldChange)
    r <- r[, .SD[1], by = ref_key]                 # best barcode per WT partner
    r[, paired_unit := u]
  })
  ref <- rbindlist(ref_list, use.names = TRUE, fill = TRUE)
  
  out <- rbind(test, ref, use.names = TRUE, fill = TRUE)
  setorder(out, log2FoldChange)
  out[]
}

# ---------------------------------------------------------------------------
# 6. plots
# ---------------------------------------------------------------------------
plot_hits <- function(gates_dt, gp_dt, selected, COMP_MAIN_db, COMP_MAIN_sg,
                      COMP_CTRL = "UT vs B", minBase = 30,
                      unit = "gene_id", folder_output) {
  
  pal <- palette_mutprof
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
                                    (grepl("ref", mutation_profile) & ref_key %in% rk)]
  
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
      scale_color_manual(values = pal) +
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
  bcs     <- unique(bcs, by = "barcode")   # a shared WT appears once per unit it controls
  setorder(bcs, ref_key, mutation_profile, stat)
  keep_cols <- intersect(c("gene_id","barcode","mutation_profile","ag_class","paired_unit",
                           "gene_name","gene","ref_key","comparison","baseMean","stat",
                           "log2FoldChange","pvalue","padj"), names(bcs))
  fwrite(bcs[, ..keep_cols], file.path(fo, "03_required_barcodes.txt"))
  
  # ---- 4 output tables ----------------------------------------------------
  # 1: full audit trail, every unit, both pools
  fwrite(cand, file.path(fo, "01_units_all.txt"))
  # 2: the decision table -- tiered hits only
  fwrite(called, file.path(fo, "02_candidates_tiered.txt"))
  # 3: 03_required_barcodes.txt (written above) -- barcode to carry forward
  # 4: barcode-level evidence behind the selection, including the WT refs that
  #    the selected units are controlled by (matched via ref_key)
  sel_rk <- .rk_of(gates, selected[[unit]], unit)
  fwrite(gates[(get(unit) %in% selected[[unit]] |
                  (grepl("ref", mutation_profile) & ref_key %in% sel_rk)) &
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
# (b)/(f) transport TAA ref_key into ref_key_3p and verify the integration at load
screen_1  <- .check_taa_integration(.transport_taa_refkey(
  fread("PEPITOPE/s7_NEW_screen_results/output_screen_1/updated_screen_results/screen1_results.txt")), "screen_1")
screen_2a <- .check_taa_integration(.transport_taa_refkey(
  fread("PEPITOPE/s7_NEW_screen_results/output_screen_2a/updated_screen_results/screen2a_results.txt")), "screen_2a")

run_pipeline(screen_dt    = screen_1,
             COMP_MAIN_db = "Cluster vs UT",
             COMP_MAIN_sg = "Single vs UT",
             pat_folder   = "P1",
             topN         = 60,
             minBase      = 100,
             gp_dt = fread("PEPITOPE/s7_NEW_screen_results/output_screen_1/updated_screen_results/gene_profiles.txt"))

run_pipeline(screen_dt    = screen_2a,
             COMP_MAIN_db = "Cluster vs UT",
             COMP_MAIN_sg = "Single vs UT",
             pat_folder   = "P5",
             topN         = 60,
             minBase      = 100,
             gp_dt = fread("PEPITOPE/s7_NEW_screen_results/output_screen_2a/updated_screen_results/gene_profiles.txt"))

run_pipeline(screen_dt    = screen_1,
             COMP_MAIN_db = "X1D3 vs UT",
             COMP_MAIN_sg = "X1D3 vs UT",
             pat_folder   = "P1_1d3_control",
             topN         = 30,
             minBase      = 100,
             gp_dt = fread("PEPITOPE/s7_NEW_screen_results/output_screen_1/updated_screen_results/gene_profiles.txt"))