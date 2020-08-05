#' Calculate response matrices using fmat multiplication
#'
prs_new <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed = 241956) {
  cmat <- get_cmat(wt)
  cmat_sqrt <- get_cmat_sqrt(wt)
  idmat <- diag(1, 3 * get_nsites(wt), 3 * get_nsites(wt))

  nsites <- get_nsites(wt)

  f <- generate_fmat(wt, nmut_per_site, mut_dl_sigma, mut_sd_min, seed)

  df2ij <- calculate_s2ij(f, idmat) %>%
    rename(df2ij = mij)
  dr2ij <- calculate_s2ij(f, cmat) %>%
    rename(dr2ij = mij)
  de2ij <- calculate_s2ij(f, cmat_sqrt) %>%
    rename(de2ij = mij)

  df2ij %>%
    inner_join(dr2ij) %>%
    inner_join(de2ij)
}

#' Calculate response matrix given fmat
#'
calculate_s2ij<- function(f, amat) {
  nsites <- dim(f)[2]
  nmut <- dim(f)[3]

  dim(f) <- c(3 * nsites, nsites * nmut)

  smat <- as.matrix(Matrix(amat) %*% Matrix(f, sparse = T))

  dim(smat) = c(3, nsites, nsites, nmut)

  s2ij <- smat^2 %>%
    apply(MARGIN = c(2, 3, 4), FUN = sum)  %>%
    apply(MARGIN = c(1, 2), FUN = mean) %>%
    matrix_to_tibble()

  s2ij
}


#' Generate matrix fmat represanting mutational scan
#'
#' f(i, j, m) is force at site i for mutation m = c(1,...nmut), at site j
#'
generate_fmat <- function(wt, nmut, mut_dl_sigma, mut_sd_min,  seed) {

  mutation = seq(nmut)
  nsites <- get_nsites(wt)
  nedges <- nrow(get_graph(wt))

  dlij <- matrix(NA, nedges, nmut)
  f <- matrix(NA, 3 * nsites, nsites * nmut)
  dim(f) <- c(3 * nsites, nsites, nmut)

  for(j in seq(nsites)) {
    for (m in mutation) {
      set.seed(seed + j * m)
      dlij[, m] <- get_delta_lij(wt, site_mut = j, mut_sd_min, mut_dl_sigma)
      f[, j, m] <- get_force(wt, dlij[, m])
    }
  }

  f

}




#' calcualte all fast response matrices
#'
prs_old <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min) {

  mutants <- get_mutants_table(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min)
  mutants %>%
    calculate_dr2ij.prs() %>%
    inner_join(calculate_df2ij.prs(mutants)) %>%
    inner_join(calculate_de2ij.prs(mutants))  %>%
    group_by(i, j) %>%
    summarise(dr2ij = mean(dr2ij),
              df2ij = mean(df2ij),
              de2ij = mean(de2ij))

}



#' calcualte stress energy change vector dvsjm
#'
#'dvsjm is dvs at site j, for mutations m = c(1, 2, 3, ..., nmut_per_site)
#'
calculate_dvsjm_new <- function(wt, nmut_per_site, mut_dl_sigma, mut_sd_min,  seed =  241956) {

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



