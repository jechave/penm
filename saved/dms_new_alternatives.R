# dms alternatives


#' Calculate double-mutational-scanning matrix
#'
dmsmat_dms.new_v2 <- function(wt, nmut, mut_dl_sigma, mut_sd_min, seed) {
  # I tested this on a protein with 100 sites and 100 mutations, it's not better than the other option, and it will take more memory...

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

  dim(dri) <- c(3 * nsites, nsites * nmut)
  dim(drj) <- c(3 * nsites, nsites * nmut)

  tic()
  dridrj <- t(dri) %*% drj
  t <- toc(quiet = T)
  t_dridrj <- t$toc - t$tic


  dim(dridrj) <- c(nsites,  nmut, nsites, nmut)

  dmsmat <- -abs(dridrj) %>%
    apply(MARGIN = c(1, 3), FUN = min)

  lst(dmsmat, t_fmat, t_dridrj)
}





