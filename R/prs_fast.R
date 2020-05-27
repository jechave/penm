#' Calcualte all fast response matrices and profiles
#'
prs.fast <- function(wt, mut_dl_sigma, mut_sd_min) {

  enm_param <- get_enm_param(wt)

  mut_param <- lst(mut_dl_sigma, mut_sd_min)

  # ij response matrices
  dfij <- calculate_dfij.fast(wt, mut_dl_sigma, mut_sd_min)

  # influence profiles
  dfj <- dfij %>%
    group_by(j) %>%
    summarise(dr2j = sum(dr2ij),
              df2j = sum(df2ij),
              de2j = 1/2 * sum(de2ij),
              dvsj = 1/2 * sum(dvsij),
              dvmj = 1/2 * sum(dvmij))

  # response profiles
  dfi <- dfij %>%
    group_by(i) %>%
    summarise(dr2i = mean(dr2ij),
              df2i = mean(df2ij),
              de2i = mean(de2ij),
              dvsi = mean(dvsij),
              dvmi = mean(dvmij))

  # mode-site response matrices
  dfnj <- calculate_dfnj.fast(wt, mut_dl_sigma, mut_sd_min)

  # mode response profiles
  dfn <- dfnj %>%
    group_by(n) %>%
    summarise(dr2n = mean(dr2nj),
              df2n = mean(df2nj),
              de2n = mean(de2nj))


  lst(enm_param, mut_param, dfi, dfj,  dfij,  dfnj)
}



# Influence profiles ---------------------------------------------------------


#' Calculate several influence profiles
#'
calculate_dfj.fast <- function(wt, mut_dl_sigma, mut_sd_min,...) {

  calculate_dr2j.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    inner_join(calculate_df2j.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_de2j.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_dvsj.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_dvmj.fast(wt, mut_dl_sigma, mut_sd_min))
}


calculate_dr2j.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_dr2ij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    group_by(j) %>%
    summarise(dr2j = sum(dr2ij))
}

calculate_df2j.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_df2ij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    group_by(j) %>%
    summarise(df2j = sum(df2ij))
}

calculate_de2j.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_de2ij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    group_by(j) %>%
    summarise(de2j = 1/2 * sum(de2ij))
}

calculate_dvsj.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_dvsij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    group_by(j) %>%
    summarise(dvsj = 1/2 * sum(dvsij))
}

calculate_dvmj.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_dvmij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    group_by(j) %>%
    summarise(dvmj = 1/2 * sum(dvmij))
}





# Response profiles -------------------------------------------------------

calculate_dfi.fast <- function(wt, mut_dl_sigma, mut_sd_min,...) {
  calculate_dr2i.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    inner_join(calculate_df2i.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_de2i.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_dvsi.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_dvmi.fast(wt, mut_dl_sigma, mut_sd_min))
}


calculate_dr2i.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_dr2ij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    group_by(i) %>%
    summarise(dr2i = mean(dr2ij))
}

calculate_df2i.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_df2ij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    group_by(i) %>%
    summarise(df2i = mean(df2ij))
}

calculate_de2i.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_de2ij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    group_by(i) %>%
    summarise(de2i =  mean(de2ij))
}

calculate_dvsi.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_dvsij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    group_by(i) %>%
    summarise(dvsi = mean(dvsij))
}

calculate_dvmi.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_dvmij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    group_by(i) %>%
    summarise(dvmi = mean(dvmij))
}



# Response matrices -------------------------------------------------------


#' Calculate mutational response matrices, site-by-site
calculate_dfij.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_dr2ij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    inner_join(calculate_df2ij.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_de2ij.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_dvsij.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_dvmij.fast(wt, mut_dl_sigma, mut_sd_min))
}


calculate_dr2ij.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  cmat <- get_cmat(wt)
  calculate_Rij.fast(wt, cmat, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)
}

calculate_df2ij.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  identity_matrix = diag(nrow = get_nsites(wt) * 3, ncol = get_nsites(wt) * 3)
  calculate_Rij.fast(wt, identity_matrix, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(df2ij = mij)
}

calculate_de2ij.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  umat <- get_umat(wt)
  evalue <- get_evalue(wt)
  cmat_sqrt <- umat %*% (sqrt(1 / evalue) * t(umat))
  calculate_Rij.fast(wt, cmat_sqrt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(de2ij = mij)
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

  dvsij %>%
    matrix_to_tibble() %>%
    rename(dvsij = mij)
}

calculate_dvmij.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  calculate_dvsij.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    inner_join(calculate_de2ij.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    mutate(dvmij = dvsij - de2ij) %>%
    select(i, j, dvmij)
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
  calculate_dr2nj.fast(wt, mut_dl_sigma, mut_sd_min) %>%
    inner_join(calculate_df2nj.fast(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_de2nj.fast(wt, mut_dl_sigma, mut_sd_min))
}


calculate_dr2nj.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector <- 1 / get_evalue(wt)
  calculate_Rnj.fast(wt, avector, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, dr2nj = mij)
}

calculate_df2nj.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector = rep(1, length(get_mode(wt)))
  calculate_Rnj.fast(wt, avector, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, df2nj = mij)
}

calculate_de2nj.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector <- sqrt(1 / get_evalue(wt))
  calculate_Rnj.fast(wt, avector, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, de2nj = mij)
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
