#' Calculate Bahar's PRS matrix analytically
#'
aprs <- function(cmat, normalise = F) {
  nsites <- nrow(cmat) / 3
  prsmat <- matrix(0, nrow = nsites, ncol = nsites)
  for (i in 1:nsites) {
    for (j in 1:nsites) {
      oi <- 3 * (i - 1) + 1
      ooi <- oi + 2
      oj <- 3 * (j - 1) + 1
      ooj <- oj + 2
      cij <- cmat[oi:ooi, oj:ooj]
      prsmat[i,j] = sum(diag(cij %*% cij)) / 3
    }
  }
  if (normalise) {
    prsmat <- diag(1 / diag(prsmat)) %*% prsmat
  }

  dfij <- t(prsmat) %>%
    matrix_to_tibble(value_name = "dr2ij")

  dfij
}
