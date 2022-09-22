#' Calculate Bahar's PRS matrix analytically
#'
#' @param cmat is a covariance matrix
#' @param normalize is a logical valuable of whether to normalize the PRS matrix (default is F)
#'
#' @returns
#'
#' @export
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

  prsmat

}
