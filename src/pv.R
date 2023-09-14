two_sided_p <- function(p) 2 * pmin(p, 1 - p)

ecdf_pvals <- function(obs, reference, thr = 0.01) {
  if (is.data.frame(reference)) {
    reference <- reference %>%
      map(ecdf)
  }
  n <- environment(reference[[1]])$nobs

  p_df <- imap(
    reference,
    ~ .x(obs[[.y]])
  ) %>%
    map(two_sided_p) %>%
    map(pmax, 1 / n) %>%
    bind_cols()

  colSums(
    apply(p_df, 1, p.adjust, "holm") < thr
  )
}
