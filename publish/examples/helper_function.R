# helper_function.R
# Helper functions for the STADyUM human CD4+ tutorial.
# Sourced by human_cd4.Rmd.

# ── TSS preparation ────────────────────────────────────────────────────────────

# Retain only the most upstream TSS per gene (strand-aware).
keep_upstream_tss <- function(tsn) {
  as.data.frame(tsn) %>%
    group_by(ensembl_gene_id, strand) %>%
    slice_min(order_by = ifelse(strand == "+", start, -start),
              n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    as_granges()
}


# ── Read-count region construction ────────────────────────────────────────────

# Build pause windows (width = kmax from TSS) and trimmed gene body regions.
# Returns list(pause, gene_body).
build_readcount_regions <- function(bw_tsn,
                                    transcripts,
                                    tsn_cutoff    = 5,
                                    gb_min_length = 6000,
                                    trim_len      = 2000,
                                    kmax          = 200) {
  bw_tsn   <- bw_tsn[bw_tsn$score >= tsn_cutoff]
  bw_pause <- promoters(bw_tsn, upstream = 0, downstream = kmax)
  rm_cols  <- intersect(c("score", "type"), colnames(mcols(bw_pause)))
  if (length(rm_cols) > 0) mcols(bw_pause)[rm_cols] <- NULL

  gngrng <- transcripts %>%
    plyranges::group_by(ensembl_gene_id) %>%
    plyranges::reduce_ranges_directed() %>%
    sort()
  bw_tts <- gngrng %>% plyranges::anchor_3p() %>% mutate(width = 1)

  idx_tts  <- match(bw_pause$ensembl_gene_id, bw_tts$ensembl_gene_id)
  idx_tsn  <- match(bw_pause$ensembl_gene_id, bw_tsn$ensembl_gene_id)
  keep     <- !is.na(idx_tts) & !is.na(idx_tsn)
  bw_tts_f <- bw_tts[idx_tts[keep]]
  bw_tsn_f <- bw_tsn[idx_tsn[keep]]
  gb       <- punion(bw_tsn_f, bw_tts_f, fill.gap = TRUE)
  gb$gene_id <- bw_tsn_f$ensembl_gene_id
  gb_filt  <- (gb[width(gb) > gb_min_length]) - trim_len

  match_idx <- match(gb_filt$gene_id, bw_pause$ensembl_gene_id)
  keep2     <- !is.na(match_idx)
  gb_filt   <- gb_filt[keep2]
  bw_pause_f <- bw_pause[match_idx[keep2]]
  names(mcols(bw_pause_f))[
    names(mcols(bw_pause_f)) == "ensembl_gene_id"] <- "gene_id"

  message("Genes retained after gene body length filter: ", length(gb_filt))
  list(pause = bw_pause_f, gene_body = gb_filt)
}


# ── Rate post-processing ───────────────────────────────────────────────────────

# Find the valley between Sharp and Broad peaks in the fkSD histogram.
find_valley_threshold <- function(x, from = 20, to = 50, bins = 100, smooth_k = 5) {
  h     <- hist(x, breaks = bins, plot = FALSE)
  cnts  <- as.numeric(stats::filter(h$counts, rep(1 / smooth_k, smooth_k), sides = 2))
  cnts[is.na(cnts)] <- h$counts[is.na(cnts)]
  idx   <- which(h$mids > from & h$mids < to)
  mins  <- findpeaks(-cnts[idx])[, 2]
  h$mids[idx[mins[which.min(cnts[idx[mins]])]]]
}

# Add quintile group labels, fkSD, and Sharp/Broad classification to a rates df.
process_sample_rate <- function(df, from = 20, to = 50, bins = 100) {
  df$betaGroup <- cut(
    df$betaAdp,
    breaks = quantile(df$betaAdp, probs = seq(0, 1, 0.2), na.rm = TRUE),
    labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
    include.lowest = TRUE
  )
  df$fkMeanGroup <- cut(
    df$fkMean,
    breaks = quantile(df$fkMean, probs = seq(0, 1, 0.2), na.rm = TRUE),
    labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
    include.lowest = TRUE
  )
  df$chiGroup <- cut(
    df$chi,
    breaks = quantile(df$chi, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE
  )
  df$fkSD    <- sqrt(df$fkVar)
  threshold  <- find_valley_threshold(df$fkSD, from, to, bins)
  message("Sharp/Broad fkSD threshold: ", round(threshold, 2))
  df$sdGroup <- ifelse(df$fkSD <= threshold, "Sharp", "Broad")
  attr(df, "threshold") <- threshold
  df
}


# ── Scatter / density plots ────────────────────────────────────────────────────

# log10(beta) vs log10(chi) point-density scatter with Spearman rho.
plot_beta_vs_chi <- function(rate_df, label_x, label_y,
                              xlim = NULL, ylim = NULL, cor_size = 5) {
  rho <- cor(log10(rate_df$betaAdp), log10(rate_df$chi),
             method = "spearman", use = "complete.obs")
  p <- ggplot(rate_df, aes(x = log10(betaAdp), y = log10(chi))) +
    geom_pointdensity(adjust = 0.3, size = 0.6) +
    scale_color_viridis_c() +
    theme_cowplot() +
    labs(x = expression(log[10] ~ beta),
         y = expression(log[10] ~ chi),
         color = "Density") +
    annotate("text", x = label_x, y = label_y,
             label = paste0("rho = ", round(rho, 2)),
             hjust = 0, size = cor_size) +
    theme(axis.title      = element_text(size = 14),
          axis.text       = element_text(size = 14),
          legend.title    = element_text(size = 10),
          legend.text     = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
  if (!is.null(xlim) || !is.null(ylim))
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  p
}

# 2D density of log10(beta) vs fkSD (sigma), dashed line at Sharp/Broad boundary.
plot_beta_vs_sigma_density <- function(rate_df,
                                        xlim      = c(-5, -2),
                                        ylim      = c(0, 70),
                                        bins      = 6,
                                        alpha     = 0.8,
                                        text_size = 14,
                                        cor_size  = 5) {
  label_x <- xlim[1] + 0.01 * diff(xlim)
  label_y <- ylim[2] - 0.01 * diff(ylim)

  sd_threshold <- NULL
  if ("sdGroup" %in% names(rate_df)) {
    sharp_max <- max(rate_df$fkSD[rate_df$sdGroup == "Sharp"], na.rm = TRUE)
    broad_min <- min(rate_df$fkSD[rate_df$sdGroup == "Broad"], na.rm = TRUE)
    if (is.finite(sharp_max) && is.finite(broad_min))
      sd_threshold <- (sharp_max + broad_min) / 2
  }

  rho <- cor(log10(rate_df$betaAdp), rate_df$fkSD,
             method = "spearman", use = "complete.obs")

  p <- ggplot(rate_df, aes(x = log10(betaAdp), y = fkSD)) +
    stat_density_2d(aes(fill = after_stat(level)),
                    geom = "polygon", bins = bins, alpha = alpha) +
    annotate("text", x = label_x, y = label_y,
             label = paste0("rho = ", round(rho, 2)),
             hjust = 0, size = cor_size) +
    scale_fill_viridis_c(name = "Density") +
    labs(x = expression(log[10] ~ beta), y = expression(sigma)) +
    cowplot::theme_cowplot() +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme(text = element_text(size = text_size))
  if (!is.null(sd_threshold))
    p <- p + geom_hline(yintercept = sd_threshold, linetype = "dashed")
  p
}

# fkMean (mu) vs fkSD (sigma) point-density scatter with Spearman rho.
plot_fkMean_vs_fkSD <- function(rate_df, cor_size = 5) {
  rho <- cor(rate_df$fkMean, rate_df$fkSD,
             method = "spearman", use = "complete.obs")
  ggplot(rate_df, aes(x = fkMean, y = fkSD)) +
    geom_pointdensity(adjust = 0.3, size = 0.6) +
    scale_color_viridis_c() +
    theme_cowplot() +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    labs(x = expression(mu), y = expression(sigma), color = "Density") +
    annotate("text", x = -Inf, y = Inf,
             label = paste0("rho = ", round(rho, 2)),
             hjust = -0.1, vjust = 1.2, size = cor_size) +
    coord_cartesian(clip = "off") +
    theme(axis.title      = element_text(size = 14),
          axis.text       = element_text(size = 14),
          legend.title    = element_text(size = 10),
          legend.text     = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"),
          plot.margin     = margin(t = 8, r = 5, b = 5, l = 5))
}


# ── Nucleosome signal extraction ───────────────────────────────────────────────

# Parse a Micro-C 1D contact BED file into a GRanges with seqinfo attached.
get_1d_signal_from_contact_df <- function(contacts_df_path, seq_info) {
  contacts_df <- read.table(
    contacts_df_path, sep = "\t", header = FALSE,
    col.names = c("chrom", "start", "end", "name", "score", "strand", "tss"),
    stringsAsFactors = FALSE
  )
  gr <- makeGRangesFromDataFrame(
    contacts_df,
    seqnames.field = "chrom", start.field = "start", end.field = "end",
    strand.field = "strand", keep.extra.columns = TRUE,
    starts.in.df.are.0based = TRUE
  )
  seqinfo(gr) <- seq_info[seqnames(seqinfo(gr))]
  gr
}

loess_smooth <- function(y) {
  x_pos <- seq(from = -1000, to = 1000, length.out = 2000)
  predict(loess(y ~ x_pos, span = 0.05), x_pos)
}

# Extract 2000 bp TSS-centered windows matched to genes in a STADyUM object.
get_tid_from_stadyum_obj <- function(stadyum_obj) {
  tsn <- stadyum_obj@pauseRegions %>% anchor_5p() %>% mutate(width = 1)
  tid <- resize(tsn, width = 2000, fix = "center")
  rate_df     <- rates(stadyum_obj)
  matched_idx <- match(rate_df$geneId, tid$gene_id)
  tid[matched_idx[!is.na(matched_idx)], ]
}

# Build a smoothed nucleosome signal matrix from a STADyUM object and 1D signal file.
# Returns list(raw, smooth) — one row per gene, one column per bp position.
get_matrix_from_nucleosome_signal <- function(rate_object, signal_file,
                                               seq_info,
                                               smooth_func = loess_smooth) {
  gr           <- get_1d_signal_from_contact_df(signal_file, seq_info)
  tid          <- get_tid_from_stadyum_obj(rate_object)
  score_matrix <- ScoreMatrix(gr, windows = tid,
                               weight.col = "score", strand.aware = TRUE)
  score_matrix <- as.matrix(score_matrix@.Data)
  list(raw    = score_matrix,
       smooth = t(apply(score_matrix, 1, smooth_func)))
}


# ── Nucleosome profile plotting ────────────────────────────────────────────────

# Line profiles faceted by class_col, colored by group_col.
# count_matrix rows must match rows of rate_df.
plot_beta_group_line_faceted <- function(count_matrix,
                                          rate_df,
                                          group_col    = "betaGroup",
                                          class_col    = "chiGroup",
                                          range        = NULL,
                                          colors       = NULL,
                                          x_label      = c(0, 200),
                                          shift        = 1001,
                                          legend_title = expression(beta ~ group),
                                          y_label      = "Mean Micro-C signal") {
  df        <- as.data.frame(count_matrix)
  df$group  <- rate_df[[group_col]]
  df$class  <- factor(rate_df[[class_col]])

  grp_means <- df %>%
    group_by(class, group) %>%
    summarise(across(starts_with("V"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  df_long <- grp_means %>%
    pivot_longer(starts_with("V"), names_to = "pos_str", values_to = "value") %>%
    mutate(position = as.numeric(str_remove(pos_str, "V")) - shift)

  if (!is.null(range))
    df_long <- df_long %>% filter(position >= range[1], position <= range[2])
  df_long <- df_long %>% filter(!class %in% c("Unknown", "Other") & !is.na(class))

  if (is.null(colors))
    colors <- RColorBrewer::brewer.pal(6, "Blues")[2:6]

  ggplot(df_long, aes(x = position, y = value, color = group)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~class, ncol = 6) +
    labs(x = "Distance from TSS (bp)", y = y_label, color = legend_title) +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = x_label) +
    scale_y_continuous(n.breaks = 4) +
    theme_cowplot() +
    theme(panel.grid       = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black",
                                          linewidth = 0.4),
          legend.position  = "right",
          axis.ticks       = element_line(color = "black", linewidth = 0.3),
          axis.line        = element_line(color = "black", linewidth = 0.3),
          axis.text        = element_text(size = 11),
          strip.text       = element_text(size = 9))
}

# Wrapper applying the standard 0–1000 bp range and 0/500 axis-tick defaults.
make_microc_profile_plot <- function(ns_matrix, group_col, class_col, ...) {
  plot_beta_group_line_faceted(
    ns_matrix$smooth,
    ns_matrix$rate_df,
    range     = c(0, 1000),
    x_label   = c(0, 500),
    group_col = group_col,
    class_col = class_col,
    ...
  )
}


# ── Histone modification helpers ───────────────────────────────────────────────

# Build a 500 bp promoter GRanges centered on each TSS, joined with rate info.
get_promoter_df_from_stadyum_obj <- function(stadyum_obj, promoter_len = 500) {
  tsn <- stadyum_obj@pauseRegions %>% anchor_5p() %>% mutate(width = 1)
  tsn <- tsn[tsn$gene_id %in% rates(stadyum_obj)$geneId]
  promoter_grng <- resize(tsn, promoter_len, fix = "center")
  promoter_df   <- as.data.frame(promoter_grng) %>%
    left_join(rates(stadyum_obj), by = c("gene_id" = "geneId"))
  promoter_df$pauseCount <- sapply(promoter_df$actualPauseSiteCounts, sum)
  promoter_df <- promoter_df %>% dplyr::select(-actualPauseSiteCounts)
  list(df = promoter_df, grang = promoter_grng)
}

# Overlap promoter windows with each bigWig in meta_sel; add mean signal columns.
get_histone_modification_sum <- function(promoter_df, promoter_grng, meta_sel) {
  for (i in seq_len(nrow(meta_sel))) {
    bw_gr <- import.bw(meta_sel$bw[i])
    hits  <- findOverlaps(promoter_grng, bw_gr)
    scores <- data.frame(
      idx   = queryHits(hits),
      score = mcols(bw_gr)$score[subjectHits(hits)]
    ) %>%
      group_by(idx) %>%
      summarize(mean_score = mean(score), .groups = "drop")
    promoter_df[, meta_sel$Target[i]] <- scores$mean_score
  }
  promoter_df
}

my_boxplot_style <- function() {
  list(
    geom_boxplot(outlier.shape = NA, linewidth = 0.3),
    scale_fill_brewer(palette = "Blues"),
    theme_cowplot(),
    theme(legend.position  = "none",
          strip.background = element_rect(fill = "white", colour = "black",
                                          linewidth = 0.4),
          axis.title.y     = element_text(size = 12),
          axis.title.x     = element_text(size = 12),
          axis.text        = element_text(size = 11),
          strip.text       = element_text(size = 9),
          axis.line        = element_line(linewidth = 0.3),
          axis.ticks       = element_line(linewidth = 0.3))
  )
}

# Boxplot of log2(histone signal) by beta quintile, faceted by chi tertile.
plot_histone_boxplot <- function(df, histone, ylim, label_x) {
  y_labels <- list(
    H3K4me3 = expression(log[2] ~ "(H3K4me3)"),
    H3K27ac = expression(log[2] ~ "(H3K27ac)")
  )
  ggplot(df, aes(x = betaGroup, y = log2(.data[[histone]]), fill = betaGroup)) +
    facet_wrap(~chiGroup, nrow = 1, ncol = 3) +
    stat_compare_means(label.x = label_x,
                       mapping = aes(label = gsub("<", "p < ",
                                                  after_stat(p.format)))) +
    labs(x    = expression(beta ~ "Groups"),
         y    = y_labels[[histone]],
         fill = expression(beta ~ group)) +
    coord_cartesian(ylim = ylim) +
    my_boxplot_style() +
    theme(axis.title.x    = element_blank(),
          axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          legend.position = "right")
}
