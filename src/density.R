simulate_from_density <- function(vec, N = 1e5) {
  rnorm(N, sample(vec, size = N, replace = TRUE), bw.SJ(vec))
}
