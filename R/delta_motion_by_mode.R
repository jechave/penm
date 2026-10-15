#' Compare two protein ensembles mode by mode (mode-dependent profiles)
#'
#' Given two proteins, compare their conformational ensembles (fluctuation patterns),
#' the proteins'  \code{cmat}, and normal modes (Principal components) are assumed known.
#'
#' (This version works only for wt and mut with no indels)
#'
#' @param wt A protein, with normal modes: they are the basis of the profiles
#' @param mut A second protein, with normal modes (a genm mutant gets them from
#'   [set_enm_nma()])
#' @param beta inverse temperature \eqn{1 / RT}{1 / RT}, in mol/kcal (see
#'   [beta_boltzmann()]); taken by `delta_motion_dmsfn` only, whose result is in Å²
#'
#' @return A vector \code{(x_n)} of size \code{nmodes}, where \code{x_n} is the property compared, for mode n.
#'
#' @seealso [delta_motion_by_site] for the same comparisons resolved by site rather
#'   than by mode. [delta_structure_by_mode] compares structure rather than motion.
#'   [get_mutant_site()] produces the `mut` argument; [set_enm()] the `wt`.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, ensemble = 7)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1)
#'
#' dmsfn <- delta_motion_dmsfn(wt, mut)
#' length(dmsfn)                       # one value per mode
#'
#' # All zero here: "lfenm" leaves kmat, and so the normal modes, unchanged.
#' range(dmsfn)
#'
#' # Under "genm" the mutation changes the network, and so the motion.
#' # A genm mutant needs its normal modes added first.
#' wt_genm <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
#'                    mut_model = "genm", d_max_graph = 14)
#' mut_genm <- set_enm_nma(get_mutant_site(wt_genm, site_mut = 11, mutation = 1))
#' range(delta_motion_dmsfn(wt_genm, mut_genm))
#'
#' @name delta_motion_by_mode
#'
NULL


#' @rdname delta_motion_by_mode
#'
#' @details `delta_motion_dmsfn` returns mode-dependent profile of changes of mean-square fluctuations \eqn{\delta \sigma_n^2}{delta sigma_n^2}.
#' This version calculates fluctuations along normal modes of wt, it's independent of possible assignment issues.
#'
#'
#' @export
#'
delta_motion_dmsfn <- function(wt, mut, beta = beta_boltzmann()) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  msf_wt <- get_msf_mode(wt, beta)
  msf_mut <- diag(t(get_umat(wt)) %*% (get_cmat(mut, beta) %*% get_umat(wt)))
  dmsf = msf_mut - msf_wt
  dmsf
}

#' @rdname delta_motion_by_mode
#'
#' @details `delta_motion_dhn` returns mode-dependent profile of entropy differences \eqn{\delta H_n}{delta H_n}
#' This version calculates fluctuations along normal modes of wt, it's independent of possible assignment issues.
#'
#' @export
#'
delta_motion_dhn <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  msf_wt <- get_msf_mode(wt)
  msf_mut <- diag(t(get_umat(wt)) %*% (get_cmat(mut) %*% get_umat(wt)))
  ha <- 1/2 * log(2 * pi * exp(1) * msf_wt)
  hb <- 1/2 * log(2 * pi * exp(1) * msf_mut)
  dhn <- hb - ha
  dhn
}



#' @rdname delta_motion_by_mode
#'
#' @details `delta_motion_rwsipn` returns mode-dependent profile of rwsip similarity
#'
#' @export
#'
delta_motion_rwsipn <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  s2wt <- get_msf_mode(wt)
  s2mut <- get_msf_mode(mut)
  wmat <- tcrossprod(s2wt, s2mut)
  wmat <- wmat / tr(wmat)
  overlap <- crossprod(get_umat(wt), get_umat(mut))
  wsip_mat <- wmat * overlap^2
  rwsipn <- sqrt(rowSums(wsip_mat))
  rwsipn
}




#' @rdname delta_motion_by_mode
#'
#' @details `delta_motion_nhn` returns mode-dependent profile of mode conservation  measure \eqn{nH_n}{nH_n}
#'
#' @export
#'
delta_motion_nhn <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  overlap <- crossprod(get_umat(wt), get_umat(mut))
  pmat <- overlap^2
  # p log(p) -> 0 as p -> 0, but computes as 0 * -Inf = NaN. Whether an overlap
  # comes out as exactly 0 rather than ~1e-40 depends on the arithmetic, so
  # leaving this unguarded made whole modes NaN, and which ones varied by
  # machine. Sum only the terms that contribute; a zero-probability mode
  # contributes zero entropy by definition, not an undefined value.
  plogp <- pmat * log(pmat)
  plogp[pmat == 0] <- 0
  hn <- -rowSums(plogp)
  nhn <- exp(hn)
  nhn
}







