#' Calculate mutation-response a mutation-response matrix, simulation methods
#'
#' Returns a mutation-response matrix
#' It uses a simulation method (averages over perturbations). For details see \doi{10.7717/peerj.11330}
#'
#' A site-by-site response matrix has elements Mij that measure the response (e.g. deformation) of site i averaged over mutations at site j.
#' A mode-by-site response matrix has elements Mnj that measure the response (e.g. deformation) along mode n averaged over mutations at site j.
#'
#'
#' It may calculate either  site-by-site or mode-by site response matrices
#' Three type of response may be calculated, "structure " (dr2ij and dr2nj), "energy" (de2ij and de2nj), and "force" (df2ij and df2nj).
#'
#'
#' @param wt is the (wild-type) protein to mutate (an object obtained using \code{set_enm})
#' @param nmut_per_site is the number of mutations per site to simulate
#' @param mut_dl_sigma is the standard deviation of a normal distribution from which edge-length perturbations are picked (LFENM model).
#' @param mut_sd_min is integer sequence-distance cutoff, only edges with \code{sdij >= mut_sd_min} are mutated
#' @param option is either "site" (default) or "mode"
#' @param response is either "structure" (default), "energy", or "force"
#'
#' @return A list containing several response matrices and the corresponding response and influence profiles
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pdb <- bio3d::read.pdb("2acy")
#' wt <- set_enm(pdb, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#' dmat <- smrs_all(wt, nmut_per_site = 10, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1024)
#' }
#'
#' @family mutscan functions
#'
smrs <- function(wt, nmut_per_site, mut_dl_sigma, mut_sd_min, option = "site", response = "structure", seed = 1024) {

  perturbations <- generate_perturbations(wt, nmut_per_site, mut_dl_sigma, mut_sd_min, seed)

  if (option == "site") {
    if (response == "structure") {
      result <- calculate_dr2ij_smrs(wt, perturbations$fmat)
    } else if (response == "energy") {
      result <- calculate_de2ij_smrs(wt, perturbations$fmat)
    } else if (response == "force") {
      result <- calculate_df2ij_smrs(wt, perturbations$fmat)
    } else if (response == "stress") {
      result <- calculate_dvsij_smrs(wt, perturbations$dlmat)
    }
    result <- matrix(pull(result, 3), get_nsites(wt), get_nsites(wt))
  } else if (option == "mode") {
    if (response == "structure") {
      result <- calculate_dr2nj_smrs(wt, perturbations$fmat)
    } else if (response == "energy") {
      result <- calculate_de2nj_smrs(wt, perturbations$fmat)
    } else if (response == "force") {
      result <- calculate_df2nj_smrs(wt, perturbations$fmat)
    }
    result <- matrix(pull(result, 3), get_nmodes(wt), get_nsites(wt))
  }
  result
}


#' Calculate mutation-response matrix using a simulation  method
#'
#' Returns several response matrices and profiles.
#'
#' It calculates two types of response matrices site-by-site and mode-by site.
#' A site-by-site response matrix has elements Mij that measure the response (e.g. deformation) of site i averaged over mutations at site j.
#' A mode-by-site response matrix has elements Mnj that measure the response (e.g. deformation) along mode n averaged over mutations at site j.
#' Response profiles (i.e. sums over columns), and influence profiles (i.e., sums over rows) are also returned.
#' It uses a simulation method (calculates responses for various instances of forces, then calculates averages)
#' Three type of response are calculated, "structure " (dr2ij and dr2nj), "energy" (de2ij and de2nj), and "force" (df2ij and df2nj).
#'
#' For details see \doi{10.7717/peerj.11330}
#'
#' @param wt is the (wild-type) protein to mutate (an object obtained using \code{set_enm})
#' @param nmut_per_site is the number of mutations per site to simulate
#' @param mut_model is the mutational model, which may be "lfenm" or "sclfenm"
#' @param mut_dl_sigma is the standard deviation of a normal distribution from which edge-length perturbations are picked (LFENM model).
#' @param mut_sd_min is integer sequence-distance cutoff, only edges with \code{sdij >= mut_sd_min} are mutated
#' @param seed is an integer seed to set_seed before generating random mutations
#'
#' @return A list containing several response matrices and the corresponding response and influence profiles
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' pdb <- bio3d::read.pdb("2acy")
#' wt <- set_enm(pdb, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#' dmat <- smrs_all(wt, nmut_per_site = 10, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1024)
#' }
#'
#' @family mutscan functions
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
#' @noRd
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
#' @noRd
#'
calculate_df2ij_smrs <- function(wt, fmat) {
  calculate_s2ij(fmat, diag(nrow(fmat))) %>%
    matrix_to_tibble() %>%
    rename(df2ij = mij)
}


#' Calculate energy response matrix de2(i, j) using method "new"
#'
#' @noRd
#'
calculate_de2ij_smrs <- function(wt, fmat) {
  calculate_s2ij(fmat, get_cmat_sqrt(wt)) %>%
    matrix_to_tibble() %>%
    rename(de2ij = mij)
}


#' Calculate structure response matrix dr2(i, j) using method "new"
#'
#' @noRd
#'
calculate_dr2ij_smrs <- function(wt, fmat) {
  calculate_s2ij(fmat, get_cmat(wt)) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)
}



#' Calculate response matrices given fmat
#'
#' @noRd
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
#' @noRd
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

#' Calculate stress energy of a site
#'
#' @noRd
#'
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
#'dvsjm is the stress energy change dvs due to mutation at site j, for mutations m = c(1, 2, 3, ..., nmut_per_site)
#'
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
