#' Plot heatmap with glmnet non-zero coefficients.
#' @param res_tab Results table with columns iteration, coef (list of tables),
#' task_id, and task_name.
#' @param feature_df A `data.frame` with features to plot side by side.
#' @param fmin include features present at least `fmin` times.
#' @param ttl Title
glmnet_coef_heatmap <- function(res_tab, feature_df, fmin = 1, pheno_col = NULL, ttl) {
  coef_tab <- res_tab %>%
    select(iteration, coef) %>%
    unnest(coef) %>%
    filter(feature != "(Intercept)") %>%
    group_by(feature) %>%
    filter(n() >= fmin) %>%
    mutate(strength = abs(median(coef))) %>%
    ungroup()

  max_eff <- max(abs(coef_tab$coef))

  col <- circlize::colorRamp2(c(-max_eff, 0, max_eff), c("blue", "white", "red"))

  feature_names <- coef_tab %>%
    arrange(class, -strength) %>%
    distinct(feature) %>%
    pull(feature)

  pretty_feature_names <- feature_names %>%
    str_replace("^[^_]+_", "") %>%
    str_replace("^HALLMARK_", "H_")

  m_feature <- feature_df[, feature_names]
  m_feature_ca <- ComplexHeatmap::columnAnnotation(
    phenotype = feature_df$IHC_Status,
    col = list(phenotype = pheno_col),
    annotation_name_side = "left"
  )
  m_feature_hm <- ComplexHeatmap::Heatmap(as.matrix(t(m_feature)),
    name = "feature value",
    clustering_distance_rows = function(x, y) {
      if (all(is.na(x)) | all(is.na(y))) {
        return(1)
      }
      1 - cor(x, y, method = "spearman")
    },

    row_labels = pretty_feature_names,
    show_column_names = FALSE,
    top_annotation = m_feature_ca,
    width = unit(30, "mm")
  )

  m2_median <- feature_df %>%
    select(c(all_of(feature_names), IHC_Status)) %>%
    pivot_longer(-IHC_Status, values_to="value", names_to='feature') %>%
    group_by(IHC_Status, feature) %>%
    summarize(value = median(value), .groups = "drop") %>%
    pivot_wider(feature, names_from=IHC_Status, values_from=value) %>%
    as_data_frame() %>%
    column_to_rownames('feature')

  m2_median <- m2_median[feature_names,]
  m2_median <- m2_median - apply(m2_median, 1, min)
  m2_median <- m2_median / apply(m2_median, 1, max)

  stopifnot(all(rownames(m2_median)==colnames(m_feature)))


  m2_ra <- ComplexHeatmap::rowAnnotation(
    standardized_median=ComplexHeatmap::anno_points(m2_median,
                                                    pch=1:30,
                                                    gp=grid::gpar(col=pheno_col[colnames(m2_median)])))

  rownames_ra <- ComplexHeatmap::rowAnnotation(
    rn = ComplexHeatmap::anno_text(pretty_feature_names, location = unit(0, "npc"), just = "left")
  )


  iteration_hms <- coef_tab %>%
    arrange(class, -strength) %>%
    pivot_wider(-strength, names_from = "feature", values_from = "coef") %>%
    group_by(class) %>%
    nest(m = -class) %>%
    mutate(m = map(m, ~ as.matrix(t(column_to_rownames(.x, "iteration"))))) %>%
    rowwise() %>%
    mutate(hm = list(
      ComplexHeatmap::Heatmap(m,
        column_title = class,
        name = "coef",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_labels = pretty_feature_names,
        width = unit(20, "mm"),
        col = col
      )
    )) %>%
    pull(hm) %>%
    reduce(`+`)

  ComplexHeatmap::draw(m_feature_hm + iteration_hms + m2_ra + rownames_ra, column_title = ttl, merge_legend = TRUE,
                       heatmap_legend_side = "left")
}

pseudo_volc_plot <- function(df, min_freq = .8) {
  df_label <- filter(df, prop > min_freq) %>%
    mutate(i = row_number())
  df %>%
    mutate(feature = if_else(prop > min_freq, feature, NA_character_)) %>%
    ggplot(aes(median_coef, prop, label = feature, col = feature)) +
    geom_point() +
    ggrepel::geom_label_repel(
      max.overlaps = Inf,
      force = 5,
      segment.curvature = -0.1,
      segment.ncp = 3,
      segment.angle = 20,
      fill = alpha(c("white"), 0.5)
    ) +
    facet_wrap(~class) +
    viridis::scale_color_viridis(discrete = TRUE, option = "H", na.value = "grey50") +
    xlab("coefficient") +
    ylab("proportion of models")
}
