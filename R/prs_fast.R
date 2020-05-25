
# influence profiles ---------------------------------------------------------

#' Calculate energy mutational response
fast_calculate_dfj <- function(wt, mut_dl_sigma, mut_sd_min,...) {
  de2ij <- fast_calculate_de2ij(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble(value_name = "de2ij")

  dvsij <- fast_calculate_dvsij(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble(value_name = "dvsij")

  dfij <- de2ij %>%
    inner_join(dvsij) %>%
    mutate(dvmij = dvsij - de2ij)

  result <- dfij %>%
    group_by(j) %>%
    summarise(de2j = 1/2 * sum(de2ij),
              dvsj = 1/2 * sum(dvsij),
              dvmj = 1/2 * sum(dvmij))

  result

}



# Response matrices


#' Calculate structural mutational response, site analysis
fast_calculate_dfij <- function(wt, mut_dl_sigma = 0.3, mut_sd_min = 2) {
  # structural differences, site analysis
  dr2ij <- fast_calculate_dr2ij(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)

  df2ij <- fast_calculate_df2ij(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(df2ij = mij)

  de2ij <- fast_calculate_de2ij(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble(value_name = "de2ij")

  dvsij <- fast_calculate_dvsij(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble(value_name = "dvsij")

  result <- dr2ij %>%
    inner_join(df2ij) %>%
    inner_join(de2ij) %>%
    inner_join(dvsij) %>%
    mutate(dvmij = dvsij - de2ij)

  result
}


fast_calculate_dr2ij <- function(wt, mut_dl_sigma, mut_sd_min) {
  cmat <- get_cmat(wt)
  fast_calculate_resp_ij(wt, cmat, mut_dl_sigma, mut_sd_min)
}

fast_calculate_df2ij <- function(wt, mut_dl_sigma, mut_sd_min) {
  identity_matrix = diag(nrow = get_nsites(wt) * 3, ncol = get_nsites(wt) * 3)
  fast_calculate_resp_ij(wt, identity_matrix, mut_dl_sigma, mut_sd_min)
}

fast_calculate_de2ij <- function(wt, mut_dl_sigma, mut_sd_min) {
  umat <- get_umat(wt)
  evalue <- get_evalue(wt)
  cmat_sqrt <- umat %*% (sqrt(1 / evalue) * t(umat))
  fast_calculate_resp_ij(wt, cmat_sqrt, mut_dl_sigma, mut_sd_min)
}

fast_calculate_dvsij <- function(wt, mut_dl_sigma, mut_sd_min) {
  g <- get_graph(wt) %>%
    filter(abs(sdij) >= mut_sd_min) %>%
    mutate(dvsij = 1/2 * kij * mut_dl_sigma^2)

  nsites <- get_nsites(wt)
  dvsij <- matrix(0, nsites, nsites)

  nedges <- nrow(g)
  for (edge in seq(nedges)) {
    i = g$i[edge]
    j = g$j[edge]
    kij = g$kij[edge]
    dvsij[i, j] = 1 / 2 * kij * mut_dl_sigma^2
  }
  dvsij <- dvsij + t(dvsij)
  diag(dvsij) = rowSums(dvsij)

  dvsij
}



#' Calculate fast response matrix
#'
fast_calculate_resp_ij <- function(wt, amat, mut_dl_sigma, mut_sd_min) {

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



# Structure response, mode -----------------------------------------------------------

#' Calculate structural mutational response, mode analysis
fast_calculate_dfnj <- function(wt, mut_dl_sigma = 0.3, mut_sd_min = 2) {
  # structural differences, site analysis
  dr2nj <- fast_calculate_dr2nj(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, dr2nj = mij)

  de2nj <- fast_calculate_de2nj(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, de2nj = mij)

  df2nj <- fast_calculate_df2nj(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, df2nj = mij)

  dr2nj %>%
    inner_join(de2nj) %>%
    inner_join(df2nj)
}


fast_calculate_dr2nj <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector <- 1 / get_evalue(wt)
  fast_calculate_resp_nj(wt, avector, mut_dl_sigma, mut_sd_min)
}

fast_calculate_df2nj <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector = rep(1, length(get_mode(wt)))
  fast_calculate_resp_nj(wt, avector, mut_dl_sigma, mut_sd_min)
}

fast_calculate_de2nj <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector <- sqrt(1 / get_evalue(wt))
  fast_calculate_resp_nj(wt, avector, mut_dl_sigma, mut_sd_min)
}


#' Calculate fast response matrix
#'
fast_calculate_resp_nj <- function(wt, avector, mut_dl_sigma, mut_sd_min) {
  g <- get_graph(wt)
  eij  <- get_eij(wt)
  umat <- get_umat(wt)

  nsites <- get_nsites(wt)
  nmodes <- length(get_mode(wt))

  stopifnot(length(avector) == nmodes)

  # response to edge perturbations
  nedges <- nrow(g)
  rnij <- matrix(NA, nmodes, nedges)

  tumat <- t(umat)
  dim(tumat) <- c(nmodes, 3, nsites)

  for (edge in seq(nedges)) {
    i = g$i[edge]
    j = g$j[edge]
    daij <- avector * as.vector((tumat[,,j] - tumat[,,i]) %*% eij[edge, ])
    rnij[, edge] <- daij^2
  }

  # response to site perturbations
  weight <-  ifelse(g$sdij >= mut_sd_min, 1, 0) * g$kij^2 * mut_dl_sigma^2
  rnij <- t(t(rnij) * weight) # note that * mutliplies by row

  rnj <- matrix(NA, nmodes, nsites)
  for (site in seq(nsites)) {
    edges_of_site <- g$i == site | g$j == site
    rnj[, site] <- rowSums(rnij[, edges_of_site])
  }

  rnj


}
