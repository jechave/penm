#' Calculate response matrices using fmat multiplication
#'
prs.new <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed) {

  f <- generate_fmat(wt, nmut_per_site, mut_dl_sigma, mut_sd_min, seed)

  calculate_df2ij.new(wt, f) %>%
    inner_join(calculate_dr2ij.new(wt, f)) %>%
    inner_join(calculate_de2ij.new(wt, f))
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




#' Calciulate mean structural response matris dr2(i, j) using method "new"
#'
calculate_dr2ij.new <- function(wt, f) {
  calculate_s2ij(f, get_cmat(wt)) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)
}

#' Calciulate mean structural response matris de2(i, j) using method "new"
#'
calculate_de2ij.new <- function(wt, f) {
  calculate_s2ij(f, get_cmat_sqrt(wt)) %>%
    matrix_to_tibble() %>%
    rename(de2ij = mij)
}

#' Calciulate mean structural response matris df2(i, j) using method "new"
#'
calculate_df2ij.new <- function(wt, f) {
  calculate_s2ij(f, diag(nrow(f))) %>%
    matrix_to_tibble() %>%
    rename(df2ij = mij)
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
    apply(MARGIN = c(1, 2), FUN = mean)

  s2ij
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




