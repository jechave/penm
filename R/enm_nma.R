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
set_enm_nma <- function(kmat, TOL = 1.e-5) {
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


