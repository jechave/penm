#' Perform Normal Mode Analysis
#'
#' Given an enm `kmat`, perform NMA
#'
#' @param kmat The K matrix to diagonalize
#' @param TOL=1.e-5 A small value, eigenvectors with eigenvalues larger than `TOL` are discarded
#'
#' @return A list with elements \code{mode}, \code{evalue}, \code{cmat}, and \code{umat}
#'
#' @examples
#'
#' enm_anm(enm, 1.e-10)
#'
#' @export
#'
#'@family enm builders
#'
enm_nma <- function(kmat, TOL = 1.e-5) {
  # Given an enm object for which kmat has already been defined, perform NMA
  # It returns a list containing mode, evalue, cmat, umat
  # Diagonalize
  eig <- eigen(kmat, symmetric = TRUE)
  evalue <- eig$values
  umat <- eig$vectors
  modes <- evalue > TOL
  evalue <- evalue[modes]
  umat  <- umat[, modes]

  nmodes <- sum(modes)
  mode <- seq.int(from = nmodes, to = 1, by = -1)
  evalue <- evalue[mode]
  umat <- umat[, mode]
  mode <- mode[mode]


  cmat <-  umat %*% ((1 / evalue) * t(umat))

  nma <- list(
    mode = mode,
    evalue = evalue,
    cmat = cmat,
    umat = umat
  )
  nma
}


#' NMA of ENM using a basis
#'
#' Before diagonalizing `kmat`, transform into a (smaller) basis. Useful for normal-mode perturbation calculations
#'
#' @param kmat  An enm K matrix
#' @param umat0 A matrix of basis vectors (e.g. normal modes of unperturbed enm)
#' @param nbasis Number of basis vectors to use (must be smaller than `ncol(umat0)`)
#' @param TOL A number near 0 to use to discard normal modes (eval > TOL are kept)
#'
#' @return A list containing `mode, evalue, cmat, umat`
#'
#' @export
#'
#' @examples
#'
#'@family enm builders
enm_nma_basis <- function(kmat, umat0, nbasis = ncol(umat0), TOL = 1.e-5) {
  # Given an enm object for which kmat has already been defined, perform NMA
  # It returns a list containing mode, evalue, cmat, umat
  stopifnot(nbasis <= ncol(umat0)) # can't use less basis vectors than available
  nbasis_max <- ncol(umat0)
  basis <- seq(from = nbasis_max, to = nbasis_max - nbasis + 1, by = -1)
  umat0 <- umat0[,basis]

  # transform to basis representation
  kmodes <- crossprod(umat0, crossprod(kmat, umat0))

  # diagonalize
  eig <- eigen(kmodes, symmetric = TRUE)

  evalue <- eig$values
  modes <- evalue > TOL
  nmodes <- sum(modes)
  mode <- seq.int(from = nmodes, to = 1, by = -1)
  evalue <- evalue[modes]
  umat <- umat0 %*% eig$vectors[, modes]
  cmat <-  umat %*% ((1 / evalue) * t(umat))

  nma <- list(
    mode = mode,
    evalue = evalue,
    cmat = cmat,
    umat = umat
  )
  nma
}
