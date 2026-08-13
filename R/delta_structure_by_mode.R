#' Compare two protein structures in nm representation
#'
#' (This version works only for wt and mut with no indels)
#'
#' @param wt A protein object with \code{xyz} and \code{enm} defined
#' @param mut A second protein object  with \code{xyz} defined
#' @return A vector with contributions of each normal mode to the given property
#'
#' @name delta_structure_by_mode
#'
NULL


#' @rdname delta_structure_by_mode
#'
#' @details `delta_structure_dr2n` calculates de square of the mode-contributions to \eqn{\delta \mathbf{r} = \mathbf{C}\mathbf{f}}
#'
#'
#' @export
#'
delta_structure_dr2n <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dr <- as.vector(get_xyz(mut) - get_xyz(wt))
  drn <- as.vector(crossprod(get_umat(wt), dr))
  stopifnot(length(wt$nma$mode) == length(drn))
  dr2n <- drn^2
  as.vector(dr2n)
}

#' @rdname delta_structure_by_mode
#'
#' @details `delta_structure_de2n` calculates de square of the mode-contributions to \eqn{\delta \mathbf{e} = \mathbf{C}^{1/2}\mathbf{f}}
#'
#'
#' @export
#'
delta_structure_de2n <- function(wt, mut) {
  get_evalue(wt) * delta_structure_dr2n(wt, mut)
}

#' @rdname delta_structure_by_mode
#'
#' @details `delta_structure_df2n` calculates de square of the mode-contributions to the force vecgtor \eqn{\mathbf{f}}
#'
#'
#' @export
#'
delta_structure_df2n <- function(wt, mut) {
  get_evalue(wt)^2 * delta_structure_dr2n(wt, mut)
}


