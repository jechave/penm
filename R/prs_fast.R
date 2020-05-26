#' Calculate perturbation-response-scanning responses
#'
prs.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  enm_param <- get_enm_param(wt)
  mut_param <- lst(mut_dl_sigma, mut_sd_min)
  dfj <- calculate_dfj.fast(wt, mut_dl_sigma, mut_sd_min)
  dfi <- calculate_dfi.fast(wt, mut_dl_sigma, mut_sd_min)
  dfij <- calculate_dfij.fast(wt, mut_dl_sigma, mut_sd_min)
  dfnj <- calculate_dfnj.fast(wt, mut_dl_sigma, mut_sd_min)

  lst(enm_param, mut_param, dfi, dfj,  dfij,  dfnj)
}

# Influence profiles ---------------------------------------------------------


#' Calculate several influence profiles
#'
calculate_dfj.fast <- function(wt, mut_dl_sigma, mut_sd_min,...) {
  tibble(j = get_site(wt),
         dr2j = calculate_dr2j.fast(wt, mut_dl_sigma, mut_sd_min),
         df2j = calculate_df2j.fast(wt, mut_dl_sigma, mut_sd_min),
         de2j = calculate_de2j.fast(wt, mut_dl_sigma, mut_sd_min),
         dvsj = calculate_dvsj.fast(wt, mut_dl_sigma, mut_sd_min),
         dvmj = calculate_dvmj.fast(wt, mut_dl_sigma, mut_sd_min))
}


calculate_dr2j.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  colSums(calculate_dr2ij.fast(wt, mut_dl_sigma, mut_sd_min))
}

calculate_df2j.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  colSums(calculate_df2ij.fast(wt, mut_dl_sigma, mut_sd_min))
}


calculate_de2j.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  0.5 * colSums(calculate_de2ij.fast(wt, mut_dl_sigma, mut_sd_min))
}

calculate_dvsj.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  0.5 * colSums(calculate_dvsij.fast(wt, mut_dl_sigma, mut_sd_min))
}

calculate_dvmj.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  0.5 * colSums(calculate_dvmij.fast(wt, mut_dl_sigma, mut_sd_min))
}


# Response profiles -------------------------------------------------------

calculate_dfi.fast <- function(wt, mut_dl_sigma, mut_sd_min,...) {
  tibble(i = get_site(wt),
         dr2i = calculate_dr2i.fast(wt, mut_dl_sigma, mut_sd_min),
         df2i = calculate_df2i.fast(wt, mut_dl_sigma, mut_sd_min),
         de2i = calculate_de2i.fast(wt, mut_dl_sigma, mut_sd_min),
         dvsi = calculate_dvsi.fast(wt, mut_dl_sigma, mut_sd_min),
         dvmi = calculate_dvmi.fast(wt, mut_dl_sigma, mut_sd_min))
}


calculate_dr2i.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  rowMeans(calculate_dr2ij.fast(wt, mut_dl_sigma, mut_sd_min))
}

calculate_df2i.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  rowMeans(calculate_df2ij.fast(wt, mut_dl_sigma, mut_sd_min))
}

calculate_de2i.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  rowMeans(calculate_de2ij.fast(wt, mut_dl_sigma, mut_sd_min))
}

calculate_dvsi.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  rowMeans(calculate_dvsij.fast(wt, mut_dl_sigma, mut_sd_min))
}

calculate_dvmi.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  rowMeans(calculate_dvmij.fast(wt, mut_dl_sigma, mut_sd_min))
}





# Response matrices -------------------------------------------------------


#' Calculate structural mutational response, site analysis
calculate_dfij.fast <- function(wt, mut_dl_sigma = 0.3, mut_sd_min = 2) {
  # structural differences, site analysis
  dr2ij <- calculate_dr2ij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)

  df2ij <- calculate_df2ij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(df2ij = mij)

  de2ij <- calculate_de2ij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble(value_name = "de2ij")

  dvsij <- calculate_dvsij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble(value_name = "dvsij")

  dvmij <- calculate_dvmij.fast(wt, mut_dl_sigma, mut_sd_min)  %>%
    matrix_to_tibble(value_name = "dvmij")

  result <- dr2ij %>%
    inner_join(df2ij) %>%
    inner_join(de2ij) %>%
    inner_join(dvsij) %>%
    inner_join(dvmij)

  result
}


calculate_dr2ij.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  cmat <- get_cmat(wt)
  calculate_Rij.fast(wt, cmat, mut_dl_sigma, mut_sd_min)
}

calculate_df2ij.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  identity_matrix = diag(nrow = get_nsites(wt) * 3, ncol = get_nsites(wt) * 3)
  calculate_Rij.fast(wt, identity_matrix, mut_dl_sigma, mut_sd_min)
}

calculate_de2ij.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  umat <- get_umat(wt)
  evalue <- get_evalue(wt)
  cmat_sqrt <- umat %*% (sqrt(1 / evalue) * t(umat))
  calculate_Rij.fast(wt, cmat_sqrt, mut_dl_sigma, mut_sd_min)
}

calculate_dvsij.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
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

calculate_dvmij.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_dvsij.fast(wt, mut_dl_sigma, mut_sd_min) - calculate_de2ij.fast(wt, mut_dl_sigma, mut_sd_min)
}



#' Calculate fast response matrix
#'
calculate_Rij.fast <- function(wt, amat, mut_dl_sigma, mut_sd_min) {

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
calculate_dfnj.fast <- function(wt, mut_dl_sigma = 0.3, mut_sd_min = 2) {
  # structural differences, site analysis
  dr2nj <- calculate_dr2nj.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, dr2nj = mij)

  de2nj <- calculate_de2nj.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, de2nj = mij)

  df2nj <- calculate_df2nj.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, df2nj = mij)

  dr2nj %>%
    inner_join(de2nj) %>%
    inner_join(df2nj)
}


calculate_dr2nj.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector <- 1 / get_evalue(wt)
  calculate_Rnj.fast(wt, avector, mut_dl_sigma, mut_sd_min)
}

calculate_df2nj.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector = rep(1, length(get_mode(wt)))
  calculate_Rnj.fast(wt, avector, mut_dl_sigma, mut_sd_min)
}

calculate_de2nj.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector <- sqrt(1 / get_evalue(wt))
  calculate_Rnj.fast(wt, avector, mut_dl_sigma, mut_sd_min)
}


#' Calculate fast response matrix
#'
calculate_Rnj.fast <- function(wt, avector, mut_dl_sigma, mut_sd_min) {
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
