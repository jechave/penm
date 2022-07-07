# Helpers for mutscan sub-package

#' Calculate \eqn{\mathbf{K}^{1/2}}
#'
#' Calculates the (matrix) square root of the network's matrix \code{kmat}
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nsites, the matrix sqrt of \code{kmat}
#'
#' @noRd
#'
#' @family mutscan helpers
get_kmat_sqrt <- function(prot) {
  evalue <- get_evalue(prot)
  umat <- get_umat(prot)
  kmat_sqrt <- umat %*% (sqrt(evalue) * t(umat))
  kmat_sqrt
}


#' Calculate \eqn{\mathbf{C}^{1/2}}
#'
#' Calculates the (matrix) square root of the covariance matrix \code{cmat}
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nsites, the matrix sqrt of \code{cmat}
#'
#' @noRd
#'
#' @family mutscan helpers
get_cmat_sqrt <- function(prot) {
  evalue <- get_evalue(prot)
  umat <- get_umat(prot)
  cmat_sqrt <- umat %*% (sqrt(1 / evalue) * t(umat))
  cmat_sqrt
}

