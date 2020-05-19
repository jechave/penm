#' Calculate structural mutational response, site analysis
fast_delta_structure_site <- function(wt, mut_dl_sigma = 0.3, mut_sd_min = 2) {
  # structural differences, site analysis
  dr2ij <- calculate_fast_dr2ij(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)

  de2ij <- calculate_fast_de2ij(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(de2ij = mij)

  df2ij <- calculate_fast_df2ij(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(df2ij = mij)

  dr2ij %>%
    inner_join(de2ij) %>%
    inner_join(df2ij)
}


calculate_fast_dr2ij <- function(wt, mut_dl_sigma, mut_sd_min) {
  cmat <- get_cmat(wt)
  calculate_fast_response_matrix(wt, cmat, mut_dl_sigma, mut_sd_min)
}

calculate_fast_df2ij <- function(wt, mut_dl_sigma, mut_sd_min) {
  identity_matrix = diag(nrow = get_nsites(wt) * 3, ncol = get_nsites(wt) * 3)
  calculate_fast_response_matrix(wt, identity_matrix, mut_dl_sigma, mut_sd_min)
}

calculate_fast_de2ij <- function(wt, mut_dl_sigma, mut_sd_min) {
  umat <- get_umat(wt)
  evalue <- get_evalue(wt)
  cmat_sqrt <- umat %*% (sqrt(1 / evalue) * t(umat))
  calculate_fast_response_matrix(wt, cmat_sqrt, mut_dl_sigma, mut_sd_min)
}




#' Calculate fast response matrix
#'
calculate_fast_response_matrix <- function(wt, amat, mut_dl_sigma, mut_sd_min) {

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
