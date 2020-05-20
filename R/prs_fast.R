
# Energy response ---------------------------------------------------------


#' Calculate energy mutational response
fast_delta_energy <- function(wt, mut_dl_sigma, mut_sd_min,...) {
  # energy differences

  result <- tibble(j = get_site(wt),
         delta_v_min =  calculate_fast_delta_v_min(wt, mut_dl_sigma, mut_sd_min),
         delta_v_stress =  calculate_fast_delta_v_stress(wt, mut_dl_sigma, mut_sd_min))
  result
}

calculate_fast_delta_v_min <- function(wt, mut_dl_sigma, mut_sd_min) {
  de2ij <- calculate_fast_de2ij(wt, mut_dl_sigma, mut_sd_min)
  dvsij <- calculate_fast_dvstress_ij(wt, mut_dl_sigma, mut_sd_min)
  dvmij <- dvsij - de2ij
  dvmj <- 1/2 * colSums(dvmij)
  dvmj
}

calculate_fast_delta_v_stress <- function(wt, mut_dl_sigma, mut_sd_min) {
  dvstress_ij <- calculate_fast_dvstress_ij(wt, mut_dl_sigma, mut_sd_min)
  dvstress_j <- 1/2 * colSums(dvstress_ij)
  dvstress_j
}

calculate_fast_dvstress_ij <- function(wt, mut_dl_sigma, mut_sd_min) {
  g <- get_graph(wt) %>%
    filter(abs(sdij) >= mut_sd_min) %>%
    mutate(dv_stress_ij = 1/2 * kij * mut_dl_sigma^2)

  nsites <- get_nsites(wt)
  dv_stress_ij <- matrix(0, nsites, nsites)

  nedges <- nrow(g)
  for (edge in seq(nedges)) {
    i = g$i[edge]
    j = g$j[edge]
    kij = g$kij[edge]
    dv_stress_ij[i, j] = 1 / 2 * kij * mut_dl_sigma^2
  }
  dv_stress_ij <- dv_stress_ij + t(dv_stress_ij)
  diag(dv_stress_ij) = rowSums(dv_stress_ij)

  dv_stress_ij
}



# Structure response, site ------------------------------------------------


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



# Structure response, mode -----------------------------------------------------------

#' Calculate structural mutational response, mode analysis
fast_delta_structure_mode <- function(wt, mut_dl_sigma = 0.3, mut_sd_min = 2) {
  # structural differences, site analysis
  dr2nj <- calculate_fast_dr2nj(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, dr2nj = mij)

  de2nj <- calculate_fast_de2nj(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, de2nj = mij)

  df2nj <- calculate_fast_df2nj(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, df2nj = mij)

  dr2nj %>%
    inner_join(de2nj) %>%
    inner_join(df2nj)
}


calculate_fast_dr2nj <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector <- 1 / get_evalue(wt)
  calculate_fast_response_matrix_mode(wt, avector, mut_dl_sigma, mut_sd_min)
}

calculate_fast_df2nj <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector = rep(1, length(get_mode(wt)))
  calculate_fast_response_matrix_mode(wt, avector, mut_dl_sigma, mut_sd_min)
}

calculate_fast_de2nj <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector <- sqrt(1 / get_evalue(wt))
  calculate_fast_response_matrix_mode(wt, avector, mut_dl_sigma, mut_sd_min)
}


#' Calculate fast response matrix
#'
calculate_fast_response_matrix_mode <- function(wt, avector, mut_dl_sigma, mut_sd_min) {
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

