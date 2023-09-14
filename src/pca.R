vec2diag <- function(v) {
  x <- diag(length(v))
  diag(x) <- v
  x
}

# this is based on a suggestion by David Jitao
correct_batches <- function(m, pca, centers, sds, batch, variance = FALSE,
                            reference_batch = NULL) {
  stopifnot(is.null(reference_batch) || reference_batch %in% batch)
  stopifnot(length(batch) == ncol(m))
  n_pc <- length(centers[[1]])
  for (b in unique(batch)) {
    ind <- batch == b
    x <- predict(pca, t(m[, ind]))
    residuals <- m[, ind] - t(x %*% t(pca$rotation))
    x[, seq_len(n_pc)] <- sweep(x[, seq_len(n_pc), drop = FALSE], 2, centers[[b]])
    if (variance) {
      x[, seq_len(n_pc)] <- sweep(x[, seq_len(n_pc), drop = FALSE], 2, sds[[b]], "/")
      if (!is.null(reference_batch)) {
        x[, seq_len(n_pc)] <- sweep(x[, seq_len(n_pc), drop = FALSE], 2, sds[[reference_batch]], "*")
      }
    }
    if (!is.null(reference_batch)) {
      x[, seq_len(n_pc)] <- sweep(x[, seq_len(n_pc), drop = FALSE], 2, centers[[reference_batch]], "+")
    }
    m[, ind] <- t(x %*% t(pca$rotation)) + residuals
  }
  m
}


# compute center for each batch
find_centers <- function(pca, batch, PC = 1) {
  stopifnot(nrow(pca$x) == length(batch))
  unique(batch) %>%
    set_names() %>%
    map(~ colMeans(pca$x[batch == .x, PC, drop = FALSE]))
}

# compute sd for each batch
find_sd <- function(pca, batch, PC = 1) {
  stopifnot(nrow(pca$x) == length(batch))
  unique(batch) %>%
    set_names() %>%
    map(~ apply(pca$x[batch == .x, PC, drop = FALSE], 2, sd))
}

rotations_hm <- function(pca) {
  ComplexHeatmap::Heatmap(pca$rotation,
    cluster_rows = F,
    cluster_columns = F,
    show_column_names = F
  )
}

plot_diff <- function(sm) {
  df <- sm %>%
    bind_rows(.id = "dataset") %>%
    pivot_longer(-dataset, names_to = "PC", names_prefix = "PC", names_transform = list(PC = as.integer))

  df <- inner_join(df, df, by = "PC") %>%
    filter(dataset.x < dataset.y) %>%
    mutate(diff = (value.x - value.y)**2) %>%
    mutate(contrast = str_c(dataset.x, "-", dataset.y))

  ggplot(df) +
    aes(PC, diff) +
    geom_point() +
    facet_wrap(~contrast)
}

plot_centers <- function(pca, batch, max_pc = 200) {
  find_centers(pca, batch, seq_len(max_pc)) %>%
    plot_diff()
}

plot_sd <- function(pca, batch, max_pc = 200) {
  find_sd(pca, batch, seq_len(max_pc)) %>%
    plot_diff()
}


plot_pca <- function(p, b, PC1 = PC1, PC2 = PC2) {
  broom::tidy(p) %>%
    filter(PC < 10) %>%
    pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>%
    mutate(group = b) %>%
    ggplot() +
    aes({{ PC1 }}, {{ PC2 }}, col = group) +
    geom_point()
}

reduce_expression_matrix <- function(m, ngenes = 500) {
  tail(m[order(apply(m, 1, sd)), ], ngenes)
}

mean_variance_comparison_plot <- function(...) {
  var_names = names(list(...))
  d <- list(...) %>%
    map(function(.x) {
      if (is.matrix(.x)) .x else exprs(.x)
    }) %>%
    map(list) %>%
    as_tibble()
  d <- d %>%
    pivot_longer(everything(), names_to = "dataset") %>%
    mutate(
      mean = map(value, rowMeans),
      variance = map(value, ~ apply(.x, 1, var))
    ) %>%
    mutate(
      feature = map(mean, ~ as_tibble_col(names(.x), "feature")),
      mean = map(mean, as_tibble_col, "mean"),
      variance = map(variance, as_tibble_col, "variance")
    ) %>%
    select(-value) %>%
    unnest(c(feature, mean, variance)) %>%
    pivot_longer(c(mean, variance), names_to = "stat") %>%
    pivot_wider(names_from = dataset, values_from = c(value))

  plt <- function(x, title) {
    rmse <- sqrt(mean((x[[var_names[1]]] - x[[var_names[2]]])^2)) %>%
      round(3)
    ggplot(x) +
      aes_string(var_names[1], var_names[2]) +
      geom_point(alpha = 0.2) +
      geom_density_2d() +
      ggpubr::stat_cor() +
      ggtitle(title) +
      annotate("text",  x=Inf, y = Inf, label = str_glue("RMSE={rmse}"), vjust=1, hjust=1)
  }

  list(
    filter(d, stat == "mean") %>%
      plt("mean"),
    filter(d, stat == "variance") %>%
      plt("variance") +
      scale_x_log10() +
      scale_y_log10()
  )
}
