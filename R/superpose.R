# Superpose a prot onto a structure -------------------------------------------


#' Superpose a protein onto a structure
#'
#' Rotates and translates a protein so that its structure best fits `fixed`, in
#' the least-squares sense, with [bio3d::fit.xyz()]. Everything else that
#' depends on the protein's orientation is rotated with it: the network matrix
#' (`kmat`), the normal modes, and the edge directions an lfenm protein keeps for
#' its mutations. What does not depend on orientation, such as the edge lengths
#' and spring constants, the eigenvalues and the energy, is unchanged.
#'
#' Comparing two proteins site by site, with `delta_structure_dr2i()` for
#' example, needs them in the same orientation. An lfenm mutant is always in the
#' orientation of the protein `set_enm()` built. A genm mutant is superposed
#' onto the protein it was made from, so along a trajectory its orientation can
#' drift; superposing the tips of independent trajectories onto the wild type
#' puts them back in one orientation.
#'
#' @param prot a `prot`, of either mutational model, with or without normal
#'   modes
#' @param fixed coordinates to superpose `prot` onto: a vector of length
#'   `3 * nsites`, as returned by [get_xyz()]
#'
#' @returns `prot`, superposed onto `fixed`
#'
#' @export
#'
#' @seealso [get_mutant_site()], whose genm mutants this is most useful for;
#'   [set_enm_nma()].
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
#'               mut_model = "genm", d_max_graph = 14)
#'
#' # a short trajectory
#' tip <- wt
#' for (site in c(10, 40, 70)) tip <- get_mutant_site(tip, site_mut = site, mutation = 1)
#'
#' # superposed onto the wild type, with its normal modes added in that orientation
#' tip <- set_enm_nma(superpose_prot(tip, get_xyz(wt)))
#' sum(delta_structure_dr2i(wt, tip))
#'
superpose_prot <- function(prot, fixed) {
  nsites <- get_nsites(prot)
  fixed <- as.vector(fixed)
  if (length(fixed) != 3 * nsites) stop("fixed must have length 3 * nsites = ", 3 * nsites)

  mobile <- get_xyz(prot)
  all_coordinates <- seq_along(fixed)
  xyz <- as.vector(bio3d::fit.xyz(fixed = fixed, mobile = mobile,
                                  fixed.inds = all_coordinates,
                                  mobile.inds = all_coordinates))

  # the rotation fit.xyz applied, as a 3 x 3 matrix, and as the 3N x 3N matrix
  # that applies it to every node
  rotation <- superposition_rotation(mobile, xyz)
  rotation_all <- kronecker(diag(nsites), rotation)

  prot$nodes$xyz <- xyz
  prot$kmat <- rotation_all %*% get_kmat(prot) %*% t(rotation_all)
  has_modes <- !identical(prot$nma, NA)
  if (has_modes) {
    prot$nma$umat <- canonical_sign(rotation_all %*% prot$nma$umat)
    prot$nma$cmat <- rotation_all %*% prot$nma$cmat %*% t(rotation_all)
  }
  # one row per edge: a row vector e rotates to e t(rotation)
  has_eij <- !is.null(prot$internal$eij)
  if (has_eij) prot$internal$eij <- prot$internal$eij %*% t(rotation)
  prot
}


#' The rotation that takes one structure onto another
#'
#' `moved` is `mobile` rotated and translated, as [bio3d::fit.xyz()] returns
#' it. Relative to their centroids, `moved = rotation %*% mobile` node by node,
#' so `rotation` is the least-squares solution of that, which is exact up to
#' rounding.
#'
#' @param mobile coordinates before the superposition, length `3 * nsites`
#' @param moved the same coordinates after it
#'
#' @returns the 3 x 3 rotation matrix
#'
#' @noRd
#'
superposition_rotation <- function(mobile, moved) {
  before <- matrix(mobile, nrow = 3)          # column k is node k
  after <- matrix(moved, nrow = 3)
  before <- before - rowMeans(before)         # relative to the centroid
  after <- after - rowMeans(after)
  after %*% t(before) %*% solve(before %*% t(before))
}
