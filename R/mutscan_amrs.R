#' Calculate mutation-response a mutation-response matrix analitically
#'
#' Returns a mutation-response matrix
#' It uses an analytical method (closed formulas). For details see \doi{10.7717/peerj.11330}
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
#' dmat <- smrs_all(wt,  mut_dl_sigma = 0.3, mut_sd_min = 1, option = "site", resonse = "structure")
#' }
#'
#' @family mutscan functions
#'
amrs <- function(wt, mut_dl_sigma, mut_sd_min, option = "site", response = "structure") {

  if (option == "site") {
    if (response == "structure") {
      result <- calculate_dr2ij_amrs(wt, mut_dl_sigma, mut_sd_min)
    } else if (response == "energy") {
      result <- calculate_de2ij_amrs(wt, mut_dl_sigma, mut_sd_min)
    } else if (response == "force") {
      result <- calculate_df2ij_amrs(wt, mut_dl_sigma, mut_sd_min)
    } else if (response == "stress") {
      result <- calculate_dvsij_amrs(wt, mut_dl_sigma, mut_sd_min)
    }
    result <- matrix(pull(result, 3), get_nsites(wt), get_nsites(wt))
  } else if (option == "mode") {
    if (response == "structure") {
      result <- calculate_dr2nj_amrs(wt, mut_dl_sigma, mut_sd_min)
    } else if (response == "energy") {
      result <- calculate_de2nj_amrs(wt, mut_dl_sigma, mut_sd_min)
    } else if (response == "force") {
      result <- calculate_df2nj_amrs(wt, mut_dl_sigma, mut_sd_min)
    }
    result <- matrix(pull(result, 3), get_nmodes(wt), get_nsites(wt))
  }
  result
}



#' Calculate mutation-response matrices and profiles using analytical methods
#'
#' Returns several response matrices and profiles.
#'
#' It calculates two types of response matrices site-by-site and mode-by site.
#' A site-by-site response matrix has elements Mij that measure the response (e.g. deformation) of site i averaged over mutations at site j.
#' A mode-by-site response matrix has elements Mnj that measure the response (e.g. deformation) along mode n averaged over mutations at site j.
#' Response profiles (i.e. sums over columns), and influence profiles (i.e., sums over rows) are also returned.
#' It uses an analytical method (closed formulas)
#' Three type of response are calculated, "structure " (dr2ij and dr2nj), "energy" (de2ij and de2nj), and "force" (df2ij and df2nj).
#'
#' For details see \doi{10.7717/peerj.11330}
#'
#' @param wt is the (wild-type) protein to mutate (an object obtained using \code{set_enm})
#' @param mut_dl_sigma is the standard deviation of a normal distribution from which edge-length perturbations are picked (LFENM model).
#' @param mut_sd_min is integer sequence-distance cutoff, only edges with \code{sdij >= mut_sd_min} are mutated
#'
#' @return A list containing several response matrices and the corresponding response and influence profiles
#'
#' @export
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
amrs_all <- function(wt, mut_dl_sigma, mut_sd_min) {

  enm_param <- get_enm_param(wt)

  mut_param <- lst(mut_dl_sigma, mut_sd_min)

  # ij response matrices
  dfij <- calculate_dr2ij_amrs(wt, mut_dl_sigma, mut_sd_min) %>%
    inner_join(calculate_df2ij_amrs(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_de2ij_amrs(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_dvsij_amrs(wt, mut_dl_sigma, mut_sd_min)) %>%
    mutate(dvmij = dvsij - de2ij)

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
  dfnj <- calculate_dr2nj_amrs(wt, mut_dl_sigma, mut_sd_min) %>%
    inner_join(calculate_df2nj_amrs(wt, mut_dl_sigma, mut_sd_min)) %>%
    inner_join(calculate_de2nj_amrs(wt, mut_dl_sigma, mut_sd_min))

  # mode response profiles
  dfn <- dfnj %>%
    group_by(n) %>%
    summarise(dr2n = mean(dr2nj),
              df2n = mean(df2nj),
              de2n = mean(de2nj))

  lst(dfij, dfi, dfj, dfnj, dfn, enm_param, mut_param)
}


# site-by-site response matrices -------------------------------------------------------


#' Calculate site-by-site force matrix df2(n, j) using method "fast"
#'
#' @noRd
#'
calculate_df2ij_amrs <- function(wt, mut_dl_sigma, mut_sd_min) {
  identity_matrix = diag(get_nsites(wt) * 3)
  calculate_Rij_amrs(wt, identity_matrix, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(df2ij = mij)
}

#' Calculate site-by-site energy response matrix de2(n, j) using method "fast"
#'
#' @noRd
#'
calculate_de2ij_amrs <- function(wt, mut_dl_sigma, mut_sd_min) {
  umat <- get_umat(wt)
  evalue <- get_evalue(wt)
  cmat_sqrt <- umat %*% (sqrt(1 / evalue) * t(umat))
  calculate_Rij_amrs(wt, cmat_sqrt, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(de2ij = mij)
}

#' Calculate site-by-site structure response matrix dr2(n, j) using method "fast"
#'
#' @noRd
#'
calculate_dr2ij_amrs <- function(wt, mut_dl_sigma, mut_sd_min) {
  cmat <- get_cmat(wt)
  calculate_Rij_amrs(wt, cmat, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)
}

#' Calculate site-by-site stress-energy response matrix dvs(n, j) using method "fast"
#'
#' @noRd
#'
calculate_dvsij_amrs <- function(wt, mut_dl_sigma, mut_sd_min) {
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



#' Calculate site-by-site "fast" response matrix, general
#'
#' @noRd
#'
calculate_Rij_amrs <- function(wt, amat, mut_dl_sigma, mut_sd_min) {

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


# mode-by-site response matrices -----------------------------------------------------------

#' Calculate mode-by-site force matrix df2(n, j) using method "fast"
#'
#' @noRd
#'
calculate_df2nj_amrs <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector = rep(1, length(get_mode(wt)))
  calculate_Rnj_amrs(wt, avector, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, df2nj = mij)
}

#' Calculate mode-by-site energy response matrix de2(n, j) using method "fast"
#'
#' @noRd
#'
calculate_de2nj_amrs <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector <- sqrt(1 / get_evalue(wt))
  calculate_Rnj_amrs(wt, avector, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, de2nj = mij)
}

#' Calculate mode-by-site structure response matrix dr2(n, j) using method "fast"
#'
#' @noRd
#'
calculate_dr2nj_amrs <- function(wt, mut_dl_sigma, mut_sd_min) {
  avector <- 1 / get_evalue(wt)
  calculate_Rnj_amrs(wt, avector, mut_dl_sigma, mut_sd_min) %>%
    matrix_to_tibble() %>%
    rename(n = i, dr2nj = mij)
}


#' Calculate fast response matrix, general
#'
#' @noRd
#'
calculate_Rnj_amrs <- function(wt, avector, mut_dl_sigma, mut_sd_min) {
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
