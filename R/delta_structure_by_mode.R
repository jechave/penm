#' Compare two protein structures in nm representation
#'
#' (This version works only for wt and mut with no indels)
#'
#' @param wt A protein object with \code{xyz} and \code{enm} defined
#' @param mut A second protein object  with \code{xyz} defined
#' @return A vector of size \code{nmodes} with the contribution of each normal mode
#'   to the given property.
#'
#' @seealso [delta_structure_by_site] for the same differences resolved by site
#'   rather than by mode — `dr2i` and `dr2n` are the same displacement in two bases
#'   and sum to the same total. [get_mutant_site()] produces the `mut` argument;
#'   [set_enm()] the `wt`.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, frustrated = FALSE)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1, ensemble = 7)
#'
#' dr2n <- delta_structure_dr2n(wt, mut)
#' length(dr2n)                        # one value per mode
#' sum(dr2n)                           # matches sum(delta_structure_dr2i(wt, mut))
#'
#' @name delta_structure_by_mode
#'
NULL


#' @rdname delta_structure_by_mode
#'
#' @details `delta_structure_dr2n` calculates the square of the mode-contributions to \eqn{\delta \mathbf{r} = \mathbf{C}\mathbf{f}}
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
#' @details `delta_structure_de2n` calculates the square of the mode-contributions to \eqn{\delta \mathbf{e} = \mathbf{C}^{1/2}\mathbf{f}}
#'
#'
#' @export
#'
delta_structure_de2n <- function(wt, mut) {
  get_evalue(wt) * delta_structure_dr2n(wt, mut)
}

#' @rdname delta_structure_by_mode
#'
#' @details `delta_structure_df2n` calculates the square of the mode-contributions to the force vector \eqn{\mathbf{f}}
#'
#'
#' @export
#'
delta_structure_df2n <- function(wt, mut) {
  get_evalue(wt)^2 * delta_structure_dr2n(wt, mut)
}


