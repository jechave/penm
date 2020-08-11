#' Calculate response matrices using fmat multiplication
#'
prs.new <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed) {

  perturbations <- generate_perturbations(wt, nmut_per_site, mut_dl_sigma, mut_sd_min, seed)
  fmat <- perturbations$fmat
  dlmat <- perturbations$dlmat

  calculate_df2ij.new(wt, fmat) %>%
    inner_join(calculate_dr2ij.new(wt, fmat)) %>%
    inner_join(calculate_de2ij.new(wt, fmat))
}

#' calcualte all fast response matrices and profiles
#'
prs_all.new <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed) {


  enm_param <- get_enm_param(wt)
  mut_param <- lst(nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min)

  perturbations <- generate_perturbations(wt, nmut_per_site, mut_dl_sigma, mut_sd_min, seed)
  fmat <- perturbations$fmat
  dlmat <- perturbations$dlmat

  dfij <- calculate_dr2ij.new(wt, fmat) %>%
    inner_join(calculate_df2ij.new(wt, fmat)) %>%
    inner_join(calculate_de2ij.new(wt, fmat)) %>%
    inner_join(calculate_dvsij.new(wt, dlmat)) %>%
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
  dfnj = NA
  dfn = NA
#   dfnj <- calculate_df2nj.prs(mutants) %>%
#     inner_join(calculate_de2nj.prs(mutants)) %>%
#     inner_join(calculate_dr2nj.prs(mutants))
#
#   dfn <- dfnj %>%
#     group_by(n) %>%
#     summarise(dr2n = mean(dr2nj),
#               df2n = mean(df2nj),
#               de2n = mean(de2nj)) %>%
#     ungroup()

  lst(dfij, dfi, dfj, dfnj, dfn, enm_param, mut_param)
}




#' Calculate matrices of perturbations (dlmat) and forces (fmat)
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
      dlmat[, j, m] <- get_delta_lij(wt, site_mut = j, mut_sd_min, mut_dl_sigma)
      fmat[, j, m] <- get_force(wt, dlmat[, j, m])
    }
  }

  lst(dlmat, fmat)

}




#' Calciulate mean structural response matris dr2(i, j) using method "new"
#'
calculate_dr2ij.new <- function(wt, fmat) {
  calculate_s2ij(fmat, get_cmat(wt)) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)
}

#' Calciulate mean structural response matris de2(i, j) using method "new"
#'
calculate_de2ij.new <- function(wt, fmat) {
  calculate_s2ij(fmat, get_cmat_sqrt(wt)) %>%
    matrix_to_tibble() %>%
    rename(de2ij = mij)
}

#' Calciulate mean structural response matris df2(i, j) using method "new"
#'
calculate_df2ij.new <- function(wt, fmat) {
  calculate_s2ij(fmat, diag(nrow(fmat))) %>%
    matrix_to_tibble() %>%
    rename(df2ij = mij)
}



#' Calculate response matrix given fmat
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



calculate_dvsij.new <- function(wt, dlmat) {
  # structural differences, site analysis
  stopifnot(dim(dlmat[2]) == get_nsites(wt))
  nsites <- dim(dlmat)[2]
  nmut <- dim(dlmat)[3]

  dvsijm <- matrix(NA, nsites, nsites * nmut)
  dim(dvsijm) <- c(nsites, nsites, nmut)
  for (m in seq(nmut) ) {
    for (j in seq(nsites)) {
      dvsijm[ , j, m] <- calculate_dvsi.noindel.new(wt, dlmat[ , j, m])
    }
  }
  dvsij <- dvsijm %>% apply(MARGIN = c(1, 2), FUN = mean) %>%
    matrix_to_tibble() %>%
    rename(dvsij = mij)

  dvsij
}

calculate_dvsi.noindel.new <- function(wt, dlij) {
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
calculate_dvsjm.new <- function(wt, nmut_per_site, mut_dl_sigma, mut_sd_min,  seed) {

  nsites <- get_nsites(wt)
  site = seq(nsites)
  mutation <- seq(nmut_per_site)

  dvsjm <- matrix(NA, nsites, nmut_per_site)

  for(j in site) {
    for (m in mutation) {
      set.seed(seed + j * m)
      dlij <- get_delta_lij(wt, j, mut_sd_min, mut_dl_sigma)
      dvsjm[j, m] <- sum(1 / 2 * get_graph(wt)$kij * dlij^2)
    }
  }

  dvsjm

}
