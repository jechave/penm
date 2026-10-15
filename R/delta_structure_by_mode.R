#' Compare two protein structures in nm representation
#'
#' (This version works only for wt and mut with no indels)
#'
#' @param wt A protein, with normal modes: they are the basis of the profiles
#' @param mut A second protein
#' @return A vector of size \code{nmodes} with the contribution of each normal mode
#'   to the given property: in Å² for `dr2n`, kcal/mol for `de2n`, and
#'   (kcal/(mol Å))² for `df2n`.
#'
#' @seealso [delta_structure_by_site] for the same differences resolved by site
#'   rather than by mode — `dr2i` and `dr2n` are the same displacement in two bases
#'   and sum to the same total. [get_mutant_site()] produces the `mut` argument;
#'   [set_enm()] the `wt`.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, ensemble = 7)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1)
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
#' @details `delta_structure_dr2n` calculates the square of the mode-contributions to the displacement \eqn{\delta \mathbf{r}}{delta r} from \code{wt} to \code{mut}
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
#' @details `delta_structure_de2n` calculates the square of the mode-contributions to \eqn{\delta \mathbf{e} = \mathbf{K}^{1/2}\delta\mathbf{r}}{delta e = K^(1/2) delta r}
#'
#'
#' @export
#'
delta_structure_de2n <- function(wt, mut) {
  get_evalue(wt) * delta_structure_dr2n(wt, mut)
}

#' @rdname delta_structure_by_mode
#'
#' @details `delta_structure_df2n` calculates the square of the mode-contributions to the force vector \eqn{\mathbf{f}}{f}
#'
#'
#' @export
#'
delta_structure_df2n <- function(wt, mut) {
  get_evalue(wt)^2 * delta_structure_dr2n(wt, mut)
}


