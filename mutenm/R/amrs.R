#' Calcualte response matrices and profiles, analytic mrs method
#'
amrs <- function(wt,  mut_dl_sigma, mut_sd_min) {

  # calculate dr2ij
  tic()
  dr2ij <- calcualte_dr2ij_amrs(wt, mut_dl_sigma, mut_sd_min)
  t <- toc(quiet = T)
  t_dr2ij <- t$toc - t$tic

  lst(dr2ij, t_dr2ij)
}




#' Calculate site-by-site structure response matrix dr2(n, j) using method "amrs"
#'

calcualte_dr2ij_amrs <- function(wt, mut_dl_sigma, mut_sd_min) {
  cmat <- get_cmat(wt)
  calculate_Rij_amrs(wt, cmat, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)
}


#' Calculate site-by-site analytic response matrix, general
#'

calculate_Rij_amrs <- function(wt, amat, mut_dl_sigma, mut_sd_min) {

  g <- get_graph(wt)
  eij  <- get_eij(wt)

  nsites <- get_nsites(wt)
  dim(amat) <- c(nrow(amat), 3, nsites)

  # response to edge perturbations
  nedges <- nrow(g)
  rxij = matrix(NA, nsites, nedges)
  for (edge in seq(nedges)) {
    i = g$i[edge]
    j = g$j[edge]
    daij <- (amat[,,j] - amat[,,i]) %*% eij[edge, ]
    daij <- my_as_xyz(daij)
    rxij[, edge] <- colSums(daij^2)
  }



  # response to site perturbations
  weight <-  ifelse(g$sdij >= mut_sd_min, 1, 0) * g$kij^2 * mut_dl_sigma^2
  rxij <- t(t(rxij) * weight) # note that * mutliplies by row

  rxj <- matrix(NA, nsites, nsites)
  for (site in seq(nsites)) {
    edges_of_site <- g$i == site | g$j == site
    rxj[, site] <- rowSums(rxij[, edges_of_site])
  }

  rxj


}
