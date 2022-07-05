
#' Calculate double-mutational-scanning matrix
#'
sdmrs.old <- function(wt, nmut, mut_dl_sigma, mut_sd_min, seed, option = "mean_min") {
  # I changed this to the version I have used for the superfast paper, the version I shared.
  # The difference is that instead of "min_min" I new use "max_max" and instead of "mean_min" I now use "mean_max"

  cmat <- get_cmat(wt)


  tic()
  fmati <- generate_fmat_dmrs(wt, nmut, mut_dl_sigma, mut_sd_min, 1 * seed)
  fmatj <- generate_fmat_dmrs(wt, nmut, mut_dl_sigma, mut_sd_min, 2 * seed)
  t <-  toc(quiet = T)
  t_fmat <- t$toc - t$tic


  nsites <- get_nsites(wt)

  # structural change due to mutations fmati
  dim(fmati) <- c(3 * nsites, nsites * nmut)
  dri <-  as.matrix(Matrix(cmat) %*% Matrix(fmati, sparse = T))
  dim(dri) <- c(3 * nsites, nsites, nmut)


  # structural change due to mutations fmatj
  dim(fmatj) <- c(3 * nsites, nsites * nmut)
  drj <-  as.matrix(Matrix(cmat) %*% Matrix(fmatj, sparse = T))
  dim(drj) <- c(3 * nsites, nsites, nmut)

  drj2 <- drj^2 %>%
    apply(MARGIN = c(2, 3), FUN = sum)
  dmrs_matrix <- matrix(NA_real_, nrow = nsites, ncol = nsites)

  # compensation matrices

  tic()

  for (i in seq(nsites)) {
    for (j in seq(nsites)) {
      mi <- dri[ , i, ]
      mj <- drj[ , j, ]
      mimj <- t(mi) %*% mj
      if (option == "mean_min")
        dmrs_matrix[i,j] <- -sqrt(mean(matrixStats::rowMaxs(mimj^2)))
      if (option == "min_min")
        dmrs_matrix[i,j] <- min(mimj)
    }
  }
  t <-  toc(quiet = T)
  t_dridrj <- t$toc - t$tic

  lst(dmrs_matrix, t_fmat, t_dridrj)
}
