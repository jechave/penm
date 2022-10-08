#' Compare two protein ensembles mode by mode (mode-dependent profiles)
#'
#' Given two proteins, compare their conformational ensembles (fluctuation patterns),
#' the proteins'  \code{cmat}, and normal modes (Principal components) are assumed known.
#'
#' (This version works only for wt and mut with no indels)
#'
#' @param wt A protein object with \code{enm} defined
#' @param mut A second protein object  with \code{enm} defined
#'
#' @return A vector \code{(x_n)} of size \code{nmodes}, where \code{x_n} is the property compared, for mode n.
#'
#' @name delta_motion_by_mode
#'
NULL


#' @rdname delta_motion_by_mode
#'
#' @details `delta_motion_dmsfn` returns mode-dependent profile of changes of mean-square fluctuations \eqn{\delta \sigma_n^2}.
#' This version calculates fluctuations along normal modes of wt, it's independent of possible assignment issues.
#'
#'
#' @export
#'
delta_motion_dmsfn <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  msf_wt <- get_msf_mode(wt)
  msf_mut <- diag(t(get_umat(wt)) %*% (get_cmat(mut) %*% get_umat(wt)))
  dmsf = msf_mut - msf_wt
  dmsf
}

#' @rdname delta_motion_by_mode
#'
#' @details `delta_motion_dhn` returns mode-dependent profile of entropy differences \eqn{\delta H_n}
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
#' @details `delta_motion_dmsfn.old` returns mode-dependent profile of changes of mean-square fluctuations \eqn{\delta \sigma_n^2}
#' This version assumes mode n of mut corresponds to mode n of wt (i.e. no reassignment issues)
#'
#' @export
#'
delta_motion_dmsfn.old <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dmsf = get_msf_mode(mut) - get_msf_mode(wt)
  dmsf
}


#' @rdname delta_motion_by_mode
#'
#' @details `delta_motion_dhn.old` returns mode-dependent profile of entropy differences \eqn{\delta H_n}
#' This version assumes mode n of mut corresponds to mode n of wt (i.e. no reassignment issues)
#'
#' @export
#'
delta_motion_dhn.old <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  ha <- 1/2 * log(2 * pi * exp(1) * get_msf_mode(wt))
  hb <- 1/2 * log(2 * pi * exp(1) * get_msf_mode(mut))
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
#' @details `delta_motion_nhn` returns mode-dependent profile of mode conservation  measure \eqn{nH_n}
#'
#' @export
#'
delta_motion_nhn <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  overlap <- crossprod(get_umat(wt), get_umat(mut))
  pmat <- overlap^2
  hn <- -rowSums(pmat * log(pmat))
  nhn <- exp(hn)
  nhn
}







