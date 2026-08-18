# Compare protein ensembles, site by site---------------------------------------

#' Compare the motions of two proteins site by site
#'
#' Given two proteins, compare their conformational ensembles (fluctuation patterns),
#' the proteins'  \code{cmat}, and normal modes (Principal components) are assumed known.
#'
#' (This version works only for wt and mut with no indels)
#'
#' @param wt A protein object with \code{enm} defined
#' @param mut A second protein object  with \code{enm} defined
#'
#' @return A vector \code{(x_i)} of size \code{nsites}, where \code{x_i} is the property compared, for site i.
#'
#' @seealso [delta_motion_by_mode] for the same comparisons resolved by normal mode
#'   rather than by site. [delta_structure_by_site] compares structure rather than
#'   motion. [get_mutant_site()] produces the `mut` argument; [set_enm()] the `wt`.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, frustrated = FALSE)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1, seed = 1024)
#'
#' dmsfi <- delta_motion_dmsfi(wt, mut)
#' length(dmsfi)                       # one value per site
#'
#' # Note these are all zero here: "lfenm" perturbs equilibrium edge lengths but
#' # leaves kmat — and hence cmat and the normal modes — untouched, so the motion
#' # is unchanged by construction. Use a model that rebuilds the network to see a
#' # non-zero profile.
#' range(dmsfi)
#'
#' @name delta_motion_by_site
#'
NULL

#' @rdname delta_motion_by_site
#'
#' @details `delta_motion_dmsfi` returns site-dependent profile of changes of mean-square fluctuations \eqn{\delta \sigma_i^2}
#'
#' @export
#'
delta_motion_dmsfi <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dmsf = get_msf_site(mut) - get_msf_site(wt)
  dmsf
}


#' @rdname delta_motion_by_site
#'
#' @details `delta_motion_dbhati` returns site-dependent profile of dbhat distances
#'
#' @export
#'
delta_motion_dbhati <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  nsites <- get_nsites(wt)
  cmat_wt <- get_cmat(wt)
  cmat_mut <- get_cmat(mut)
  dim(cmat_wt) <- c(3, nsites, 3, nsites)
  dim(cmat_mut) <- c(3, nsites, 3, nsites)
  dbhati <- rep(NA, nsites)
  for (i in seq(nsites)) {
    dbhati[i] <- dbhat(cmat_wt[, i, , i], cmat_mut[, i, , i])
  }
  dbhati
}


#' @rdname delta_motion_by_site
#'
#' @details `delta_motion_rwsipi` returns site-dependent profile of rwsip similarity
#'
#' @export
#'
delta_motion_rwsipi <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  nsites <- get_nsites(wt)
  cmat_wt <- get_cmat(wt)
  cmat_mut <- get_cmat(mut)
  dim(cmat_wt) <- c(3, nsites, 3, nsites)
  dim(cmat_mut) <- c(3, nsites, 3, nsites)
  rwsipi <- rep(NA, nsites)
  for (i in seq(nsites)) {
    rwsipi[i] <- rwsip(cmat_wt[, i, , i], cmat_mut[, i, , i])
  }
  rwsipi
}


#' @rdname delta_motion_by_site
#'
#' @details `delta_motion_dhi` returns site-dependent profile of dh similarity
#'
#' @export
#'
delta_motion_dhi <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  nsites <- get_nsites(wt)
  cmat_wt <- get_cmat(wt)
  cmat_mut <- get_cmat(mut)
  dim(cmat_wt) <- c(3, nsites, 3, nsites)
  dim(cmat_mut) <- c(3, nsites, 3, nsites)
  dhi <- rep(NA, nsites)
  for (i in seq(nsites)) {
    dhi[i] <- dh(cmat_wt[, i, , i], cmat_mut[, i, , i])
  }
  dhi
}


