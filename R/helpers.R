# Plotting, table and analysis helpers ------------------------------------------

# Consistent, compact HTML tables
nice_table <- function(df, digits = 3, ...) {
  df %>%
    knitr::kable(format = "html", digits = digits, ...) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE,
      position = "center"
    )
}

# Save a ggplot to figures/ and return it so it also prints in the book
save_plot <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  ggplot2::ggsave(file.path("figures", filename), plot,
                  width = width, height = height, dpi = dpi, bg = "white")
  invisible(plot)
}

# Save a ComplexHeatmap (or any grid-drawing function) to figures/
save_grid <- function(draw_fun, filename, width = 7, height = 8, res = 300) {
  png(file.path("figures", filename), width = width, height = height,
      units = "in", res = res)
  on.exit(dev.off())
  draw_fun()
}

# Format an edgeR topTags table into the column layout used across chapters
# (the DGEList carries SYMBOL / ENTREZID / GENENAME in y$genes, so topTags
# already returns them as columns)
tidy_toptags <- function(qlf) {
  edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table %>%
    tibble::as_tibble() %>%
    dplyr::left_join(contamination_flag_tbl, by = "SYMBOL") %>%
    dplyr::mutate(
      status = dplyr::case_when(
        FDR < fdr_cutoff & logFC > 0 ~ "Up",
        FDR < fdr_cutoff & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      status = factor(status, levels = c("Down", "NS", "Up"))
    ) %>%
    dplyr::relocate(SYMBOL, GENENAME, ENTREZID)
}

# Summarise numbers of DE genes at the standard thresholds
de_counts <- function(res, label) {
  tibble::tibble(
    Comparison = label,
    `Up (FDR < 0.05)` = sum(res$FDR < fdr_cutoff & res$logFC > 0),
    `Down (FDR < 0.05)` = sum(res$FDR < fdr_cutoff & res$logFC < 0),
    `Up (and ≥ 2-fold)` = sum(res$FDR < fdr_cutoff & res$logFC >= lfc_cutoff),
    `Down (and ≥ 2-fold)` = sum(res$FDR < fdr_cutoff & res$logFC <= -lfc_cutoff)
  )
}

# Volcano plot with the top up/down genes labelled
plot_volcano <- function(res, n_label = 12, title = NULL, lfc_cap = 10) {
  lab <- dplyr::bind_rows(
    res %>% dplyr::filter(status == "Up")   %>% dplyr::slice_min(PValue, n = n_label),
    res %>% dplyr::filter(status == "Down") %>% dplyr::slice_min(PValue, n = n_label)
  )
  # Fold changes are capped at +/- lfc_cap so that a few genes with near-zero
  # counts in one group do not stretch the axis (they are drawn on the edge).
  res <- dplyr::mutate(res, logFC = pmax(pmin(logFC, lfc_cap), -lfc_cap))
  lab <- dplyr::mutate(lab, logFC = pmax(pmin(logFC, lfc_cap), -lfc_cap))
  ggplot(res, aes(logFC, -log10(PValue), colour = status)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed",
               colour = "grey50", linewidth = 0.3) +
    ggrepel::geom_text_repel(
      data = lab, aes(label = SYMBOL), colour = "black", size = 3,
      fontface = "italic", max.overlaps = Inf, min.segment.length = 0,
      box.padding = 0.35, seed = 1
    ) +
    scale_colour_manual(values = de_colours, name = "FDR < 0.05") +
    labs(x = expression(log[2]~fold~change), y = expression(-log[10]~P), title = title) +
    theme(legend.position = "top")
}

# MA plot with a chosen set of labels
plot_ma <- function(res, label_genes, title = NULL) {
  ggplot(res, aes(logCPM, logFC, colour = status)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    ggrepel::geom_text_repel(
      data = dplyr::filter(res, SYMBOL %in% label_genes),
      aes(label = SYMBOL), colour = "black", size = 3, fontface = "italic",
      max.overlaps = Inf, min.segment.length = 0, box.padding = 0.35, seed = 1
    ) +
    scale_colour_manual(values = de_colours, name = "FDR < 0.05") +
    labs(x = "Average expression (logCPM)", y = expression(log[2]~fold~change),
         title = title) +
    theme(legend.position = "top")
}

# Row-z-scored heatmap with genes grouped into labelled blocks ------------------
# Generalises the two pathway/mechanism heatmap functions of the original book.
# Columns are split by sample group explicitly (rather than assuming the first
# half of the columns is one group), so it works for any design.
#
#   logcpm      : genes x samples log-CPM matrix (rownames = gene symbols)
#   gene_groups : data frame with columns SYMBOL and block, in plotting order
#   samples     : data frame with columns sample and group (factor)
grouped_zscore_heatmap <- function(logcpm, gene_groups, samples,
                                   block_name = "Panel",
                                   block_wrap = 22,
                                   row_fontsize = 8,
                                   cluster_within = FALSE) {
  gene_groups <- gene_groups %>%
    dplyr::filter(SYMBOL %in% rownames(logcpm)) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
  if (nrow(gene_groups) == 0) stop("No genes to plot.")

  mat <- logcpm[gene_groups$SYMBOL, samples$sample, drop = FALSE]
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0

  blocks <- factor(gene_groups$block, levels = unique(gene_groups$block))
  levels(blocks) <- stringr::str_wrap(levels(blocks), width = block_wrap)
  block_cols <- setNames(
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nlevels(blocks)),
    levels(blocks)
  )

  grp <- droplevels(samples$group)
  top <- ComplexHeatmap::HeatmapAnnotation(
    Group = grp,
    col = list(Group = group_colours[levels(grp)]),
    annotation_label = "Group",
    show_annotation_name = FALSE,
    annotation_legend_param = list(Group = list(labels = group_labels[levels(grp)]))
  )
  left <- ComplexHeatmap::rowAnnotation(
    block = blocks,
    col = list(block = block_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  ComplexHeatmap::Heatmap(
    z,
    name = "z-score",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
    top_annotation = top,
    left_annotation = left,
    row_split = blocks,
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    row_gap = grid::unit(1.2, "mm"),
    column_split = grp,
    column_title = sub(" (saline|HDM)$", "\n\\1", group_labels[levels(grp)]),
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_gap = grid::unit(1.2, "mm"),
    cluster_rows = cluster_within,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = row_fontsize, fontface = "italic"),
    border = TRUE
  )
}
