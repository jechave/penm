# Structure differences ---------------------------------------------------

#' Calculate site-dependent profiles of structure differences between two proteins
#'
#' This version works only for wt and mut with no indels
#'
#'
#' @param wt A protein object with \code{xyz} defined
#' @param mut A second protein object  with \code{xyz} defined
#' @param kmat_sqrt The matrix square root of the ENM K matrix, as returned by
#'   \code{\link{get_kmat_sqrt}}
#'
#' @return A vector \code{(x_i)} of size \code{nsites}, where \code{x_i} is the property compared, for site i.
#'
#' @seealso [delta_structure_by_mode] for the same differences resolved by normal
#'   mode rather than by site — `dr2i` and `dr2n` are the same displacement in two
#'   bases and sum to the same total. [get_mutant_site()] produces the `mut`
#'   argument; [set_enm()] the `wt`.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1, ensemble = 7)
#'
#' dr2i <- delta_structure_dr2i(wt, mut)
#' length(dr2i)                        # one value per site
#' which.max(dr2i)                     # site that moved most
#'
#' # de2i needs the square root of the wild-type K matrix
#' delta_structure_de2i(wt, mut, kmat_sqrt = get_kmat_sqrt(wt))[1:5]
#'
#' @name delta_structure_by_site
#'
NULL

#' @rdname delta_structure_by_site
#'
#' @details `delta_structure_dr2i` returns the square of structural difference vector \eqn{\mathbf{C}\mathbf{f}}
#'
#' @export
#'
delta_structure_dr2i <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dxyz <- my_as_xyz(mut$nodes$xyz - wt$nodes$xyz) # use c(3, nsites) representation of xyz
  dr2i <- colSums(dxyz^2)
  dr2i
}

#' @rdname delta_structure_by_site
#' @details `delta_structure_de2i` returns the square of deformation energy vector \eqn{\mathbf{C}^{1/2}\mathbf{f}}
#'
#' @export
#'
delta_structure_de2i <- function(wt, mut, kmat_sqrt) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dr <- as.vector(get_xyz(mut) - get_xyz(wt))
  de <- kmat_sqrt %*% dr
  de <- my_as_xyz(de)
  de2i <-  colSums(de^2)
  de2i
}


#' @rdname delta_structure_by_site
#' @details `delta_structure_df2i` returns the square of force vector \eqn{\mathbf{f}}
#'
#' @export
#'
delta_structure_df2i <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels

  kmat <- wt$kmat

  dxyz <- my_as_xyz(mut$nodes$xyz - wt$nodes$xyz) # use c(3, nsites) representation of xyz
  dr <- as.vector(dxyz)
  df <- kmat %*% dr
  df <- my_as_xyz(df)
  df2i <-  colSums(df^2)
  df2i
}
