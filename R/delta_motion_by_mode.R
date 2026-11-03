#' Compare the fluctuations of two proteins mode by mode
#'
#' Compares the fluctuations of `wt` and `mut` along normal modes: how much the
#' fluctuation along each mode changed, and whether the mode's direction is
#' kept. The two proteins must have the same sites and be superposed
#' ([superpose_prot()]): the functions stop if they are not, and warn if they
#' are only nearly so.
#'
#' @param wt A protein, with normal modes
#' @param mut A second protein with the same sites, with normal modes (a genm
#'   mutant gets them from [set_enm_nma()])
#'
#' @return A vector of size `nmodes`, one value per mode of `wt` (or of
#'   `reference`, for `delta_motion_dh_mode()`), softest first.
#'
#' @seealso [delta_motion_by_site] for comparisons site by site; [get_h_mode()];
#'   [delta_structure_by_mode] compares structure rather than fluctuations.
#'   [get_mutant_site()] produces the `mut` argument; [set_enm()] the `wt`.
#'
#' @examples
#' # Under "genm" a mutation changes the network, and so the fluctuations.
#' # A genm mutant needs its normal modes added first.
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
#'               mut_model = "genm", d_max_graph = 14)
#' mut <- set_enm_nma(get_mutant_site(wt, site_mut = 11, mutation = 1))
#'
#' dh <- delta_motion_dh_mode(wt, mut)
#' length(dh)                          # one value per mode
#'
#' nmix <- delta_motion_nmix_mode(wt, mut)
#' range(nmix)                         # 1 for a mode kept as it is
#'
#' @name delta_motion_by_mode
#'
NULL

#' @rdname delta_motion_by_mode
#'
#' @param reference a protein with normal modes and the same sites, whose modes
#'   are the directions along which `wt` and `mut` are compared; `wt` by
#'   default. Used by `delta_motion_dh_mode()` only.
#'
#' @details `delta_motion_dh_mode()` returns the change in the entropy of the
#' fluctuation along each mode of `reference`,
#' `get_h_mode(mut, reference) - get_h_mode(wt, reference)`:
#' \deqn{\Delta h_n = \frac{1}{2} \ln \frac{\mathbf{u}_n^T \mathbf{C}^{mut} \mathbf{u}_n}{\mathbf{u}_n^T \mathbf{C}^{wt} \mathbf{u}_n},}{
#'       Delta h_n = 1/2 ln(u_n' C^mut u_n / u_n' C^wt u_n),}
#' with \eqn{\mathbf{u}_n}{u_n} the reference's mode n, up to the small
#' correction for structures that differ described in [get_h_mode()]. It is the
#' log of the ratio of the standard deviations along the mode, and positive
#' where `mut` fluctuates more. It depends only on the fluctuation along the
#' mode's direction, not on how `mut`'s own modes are arranged, so it needs no
#' pairing of the modes of the two proteins.
#'
#' For a fixed reference, swapping `wt` and `mut` changes the sign. To compare
#' many proteins along one set of directions, pass the same `reference` to every
#' comparison. Summed over modes it is not the change in the protein's entropy:
#' each mode's entropy ignores the correlations between modes.
#'
#' @export
#'
delta_motion_dh_mode <- function(wt, mut, reference = wt) {
  stopifnot(get_pdb_site(wt) == get_pdb_site(mut)) # no indels
  h_wt <- get_h_mode(wt, reference)
  h_mut <- get_h_mode(mut, reference)
  h_mut - h_wt
}


#' @rdname delta_motion_by_mode
#'
#' @details `delta_motion_nmix_mode()` returns, for each mode n of `wt`, the
#' effective number of modes of `mut` over which it is spread. With
#' \eqn{O_{nm} = \mathbf{u}^{wt}_n \cdot \mathbf{u}^{mut}_m}{O_nm = u^wt_n . u^mut_m},
#' the squared overlaps \eqn{O^2_{nm}}{O_nm^2} over m are the weights with which
#' the direction of mode n is distributed over the modes of `mut`, and
#' \deqn{\mathrm{nmix}_n = \exp\left(-\sum_m O^2_{nm} \ln O^2_{nm}\right).}{
#'       nmix_n = exp(-sum_m O_nm^2 ln O_nm^2).}
#' It is 1 when mode n is also a mode of `mut`, and k when it is spread equally
#' over k of them. It depends only on the directions of the modes, not on their
#' variances. Modes with close variances mix easily, so a large value with a
#' `delta_motion_dh_mode()` near zero usually means a mode reshuffled among
#' others of similar variance, with the fluctuation along it unchanged. If the
#' two structures differ, the weights of each mode are rescaled to sum to one
#' (see [get_h_mode()] for why they fall slightly short).
#'
#' @export
#'
delta_motion_nmix_mode <- function(wt, mut) {
  stopifnot(get_pdb_site(wt) == get_pdb_site(mut)) # no indels
  check_superposed(wt, mut)

  overlap_squared <- squared_overlap(wt, mut)
  # row n, rescaled to sum to one: how mode n of wt is distributed over the modes of mut
  weight <- overlap_squared / rowSums(overlap_squared)

  # p log(p) -> 0 as p -> 0, but computes as 0 * -Inf = NaN. Whether a weight
  # comes out as exactly 0 rather than ~1e-40 depends on the arithmetic, so
  # leaving this unguarded made whole modes NaN, and which ones varied by
  # machine. Sum only the terms that contribute; a zero weight contributes zero
  # entropy by definition, not an undefined value.
  weight_log_weight <- weight * log(weight)
  weight_log_weight[weight == 0] <- 0
  entropy <- -rowSums(weight_log_weight)
  exp(entropy)
}
