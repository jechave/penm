


#' Calcualte all response matrices and profiles, "new" prs method
#'

smrs_all <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed) {

  enm_param <- get_enm_param(wt)
  mut_param <- lst(nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min)

  perturbations <- generate_perturbations(wt, nmut_per_site, mut_dl_sigma, mut_sd_min, seed)
  fmat <- perturbations$fmat
  dlmat <- perturbations$dlmat

  dfij <- calculate_dr2ij_smrs(wt, fmat) %>%
    inner_join(calculate_df2ij_smrs(wt, fmat)) %>%
    inner_join(calculate_de2ij_smrs(wt, fmat)) %>%
    inner_join(calculate_dvsij_smrs(wt, dlmat)) %>%
    mutate(dvmij = dvsij - de2ij)

  dfj <- dfij %>%
    group_by(j) %>%
    summarise(dr2j = sum(dr2ij),
              df2j = sum(df2ij),
              de2j = 1/2 * sum(de2ij),
              dvsj = 1/2 * sum(dvsij),
              dvmj = 1/2 * sum(dvmij)) %>%
    ungroup()

  dfi <- dfij %>%
    group_by(i) %>%
    summarise(dr2i = mean(dr2ij),
              df2i = mean(df2ij),
              de2i = mean(de2ij),
              dvsi = mean(dvsij),
              dvmi = mean(dvmij)) %>%
    ungroup()
#
#   # structural differences, mode analysis
  dfnj <- calculate_df2nj_smrs(wt, fmat) %>%
    inner_join(calculate_de2nj_smrs(wt, fmat)) %>%
    inner_join(calculate_dr2nj_smrs(wt, fmat))

  dfn <- dfnj %>%
    group_by(n) %>%
    summarise(dr2n = mean(dr2nj),
              df2n = mean(df2nj),
              de2n = mean(de2nj)) %>%
    ungroup()

  lst(dfij, dfi, dfj, dfnj, dfn, enm_param, mut_param)
}


# Generate perturbations (forces) -----------------------------------------



#' Generate matrices of perturbations (dlmat) and forces (fmat)
#'
#' dlmat(edge, j, m) is the perturbation dl at edge e due to mutation m at site j
#' fmat(i, j, m) is force at site i for mutation m = c(1,...nmut), at site j
#'

generate_perturbations <- function(wt, nmut, mut_dl_sigma, mut_sd_min,  seed) {

  mutation = seq(nmut)
  nsites <- get_nsites(wt)
  nedges <- nrow(get_graph(wt))

  dlmat <- matrix(NA, nedges, nsites * nmut)
  dim(dlmat) <- c(nedges, nsites, nmut)

  fmat <- matrix(NA, 3 * nsites, nsites * nmut)
  dim(fmat) <- c(3 * nsites, nsites, nmut)

  for(j in seq(nsites)) {
    for (m in mutation) {
      set.seed(seed + j * m)
      dlmat[, j, m] <- generate_delta_lij(wt, site_mut = j, mut_sd_min, mut_dl_sigma)
      fmat[, j, m] <- calculate_force(wt, dlmat[, j, m])
    }
  }


  lst(dlmat, fmat)

}


# site-by-site response matrices ------------------------------------------



#' Calculate force response matrix df2(i, j) using method "new"
#'

calculate_df2ij_smrs <- function(wt, fmat) {
  calculate_s2ij(fmat, diag(nrow(fmat))) %>%
    matrix_to_tibble() %>%
    rename(df2ij = mij)
}


#' Calculate energy response matrix de2(i, j) using method "new"
#'

calculate_de2ij_smrs <- function(wt, fmat) {
  calculate_s2ij(fmat, get_cmat_sqrt(wt)) %>%
    matrix_to_tibble() %>%
    rename(de2ij = mij)
}


#' Calculate structure response matrix dr2(i, j) using method "new"
#'

calculate_dr2ij_smrs <- function(wt, fmat) {
  calculate_s2ij(fmat, get_cmat(wt)) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)
}



#' Calculate response matrices given fmat
#'

calculate_s2ij<- function(fmat, amat) {
  nsites <- dim(fmat)[2]
  nmut <- dim(fmat)[3]

  dim(fmat) <- c(3 * nsites, nsites * nmut)

  smat <- as.matrix(Matrix(amat) %*% Matrix(fmat, sparse = T))

  dim(smat) = c(3, nsites, nsites, nmut)

  s2ij <- smat^2 %>%
    apply(MARGIN = c(2, 3, 4), FUN = sum)  %>%
    apply(MARGIN = c(1, 2), FUN = mean)

  s2ij
}


#' Calculate stress energy response matrix dvs(i, j) using method "new"
#'

calculate_dvsij_smrs <- function(wt, dlmat) {
  # structural differences, site analysis
  stopifnot(dim(dlmat[2]) == get_nsites(wt))
  nsites <- dim(dlmat)[2]
  nmut <- dim(dlmat)[3]

  dvsijm <- matrix(NA, nsites, nsites * nmut)
  dim(dvsijm) <- c(nsites, nsites, nmut)
  for (m in seq(nmut) ) {
    for (j in seq(nsites)) {
      dvsijm[ , j, m] <- calculate_dvsi.noindel_smrs(wt, dlmat[ , j, m])
    }
  }
  dvsij <- dvsijm %>% apply(MARGIN = c(1, 2), FUN = mean) %>%
    matrix_to_tibble() %>%
    rename(dvsij = mij)

  dvsij
}

calculate_dvsi.noindel_smrs <- function(wt, dlij) {
  gwt <- get_graph(wt)

  dvsij <-  1/2 * gwt$kij * dlij^2

  dvsij_non_zero <- !near(dvsij, 0)

  dvsij <- dvsij[dvsij_non_zero]
  i <- gwt$i[dvsij_non_zero]
  j <- gwt$j[dvsij_non_zero]
  sites_non_zero <- unique(c(i,j))

  dvsi <- rep(0, get_nsites(wt))

  for (e in seq_along(dvsij))  {
    dvsi[i[e]] = dvsi[i[e]] + dvsij[e]
    dvsi[j[e]] = dvsi[j[e]] + dvsij[e]
  }

  dvsi
}


#' calcualte stress energy change vector dvsjm
#'
#'dvsjm is stress energy change dvs due to mutation at site j, for mutations m = c(1, 2, 3, ..., nmut_per_site)
#'

calculate_dvsjm_smrs <- function(wt, nmut_per_site, mut_dl_sigma, mut_sd_min,  seed) {

  nsites <- get_nsites(wt)
  site = seq(nsites)
  mutation <- seq(nmut_per_site)

  dvsjm <- matrix(NA, nsites, nmut_per_site)

  for(j in site) {
    for (m in mutation) {
      set.seed(seed + j * m)
      dlij <- generate_delta_lij(wt, j, mut_sd_min, mut_dl_sigma)
      dvsjm[j, m] <- sum(1 / 2 * get_graph(wt)$kij * dlij^2)
    }
  }

  dvsjm

}



# mode-by-site response matrices ------------------------------------------



#' Calculate mode-by-site force response matrix df2(n, j) using method "new"
#'

calculate_df2nj_smrs <- function(wt, fmat) {
  # structural differences, mode analysis
  calculate_fnmat(wt, fmat)^2 %>%
    apply(MARGIN = c(1, 2), FUN = mean) %>%
    matrix_to_tibble() %>%
    rename(n = i, j = j, df2nj = mij)
}

#' Calculate fmat in normal mode representation
#'

calculate_fnmat <- function(wt, fmat) {
  nsites <- dim(fmat)[2]
  nmut <- dim(fmat)[3]
  nmodes <- length(get_evalue(wt))
  dim(fmat) <- c(3 * nsites, nsites * nmut)
  fn <-   as.matrix(Matrix(t(get_umat(wt))) %*% Matrix(fmat, sparse = T))
  dim(fn) <- c(nmodes, nsites, nmut)
  fn
}


#' Calculate mode-by-site energy response matrix de2(n, j) using method "new"
#'

calculate_de2nj_smrs <- function(wt, fmat) {
  # structural differences, mode analysis

  fnmat <- calculate_fnmat(wt, fmat)

  denjm <- sqrt(1 / get_evalue(wt)) * fnmat

  denjm^2  %>%
    apply(MARGIN = c(1, 2), FUN = mean) %>%
    matrix_to_tibble() %>%
    rename(n = i, j = j, de2nj = mij)

}

#' Calculate mode-by-site structure response matrix dr2(n, j) using method "new"
#'

calculate_dr2nj_smrs <- function(wt, fmat) {
  # structural differences, modr analysis

  fnmat <- calculate_fnmat(wt, fmat)

  drnjm <- 1 / get_evalue(wt) * fnmat

  drnjm^2  %>%
    apply(MARGIN = c(1, 2), FUN = mean) %>%
    matrix_to_tibble() %>%
    rename(n = i, j = j, dr2nj = mij)

}
