## Energy differences
#' Calculate energy differences between a mutant and wild type
#'
#' @param wt A protein object
#' @param mut A second protein
#' @param ideal A protein object whose conformation defines the reference (ideal) state.
#'   Defaults to \code{wt}.
#' @param pdb_site_active A vector of active-site residues, in pdb numbering (\code{resno})
#' @param beta Inverse temperature, \code{1/kT}
#'
#' @return A (scalar) energy difference between mutant and wild type.
#'
#' @seealso [delta_structure_by_site] and [delta_motion_by_site] for per-site
#'   profiles rather than scalars. [get_mutant_site()] produces the `mut` argument;
#'   [set_enm()] the `wt`.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, frustrated = FALSE)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1, ensemble = 7)
#'
#' ddg_dv(wt, mut)      # minimum-energy difference
#'
#' # Under "lfenm" the contact map — and so kmat and its spectrum — is unchanged,
#' # which makes the entropic term exactly zero. It is informative only for models
#' # that rebuild the network, such as "sclfenm".
#' ddg_tds(wt, mut)
#'
#' @name delta_energy
#'
NULL

#' @rdname delta_energy
#'
#' @details `ddg_dv` calculates the minimum-energy difference between \code{mut} and \code{wt}
#'
#' @export
#'
ddg_dv <- function(wt, mut)
  enm_v_min(mut) - enm_v_min(wt)

#' @rdname delta_energy
#'
#' @details `ddg_tds` calculates the entropic free energ difference between \code{mut} and \code{wt}
#'
#' @export
#'
ddg_tds <- function(wt, mut, beta = beta_boltzmann())
  enm_g_entropy(mut, beta) - enm_g_entropy(wt, beta)


#' @rdname delta_energy
#'
#' @details `delta_energy_dvs` calculates the ideal-conformation stress-energy difference between \code{mut} and \code{wt}
#'
#' @export
#'
delta_energy_dvs <- function(wt, mut, ideal = wt)
  calculate_vs(mut, ideal) - calculate_vs(wt, ideal)


## Activation energy


#' @rdname delta_energy
#'
#' @details `ddgact_dv` calculates the energy contribution to the change in activation energy between \code{mut} and \code{wt}
#'
#' @export
#'
ddgact_dv <- function(wt, mut, ideal = wt, pdb_site_active = NA) {
  result <- dgact_dv(mut, ideal, pdb_site_active) - dgact_dv(wt, ideal, pdb_site_active)
  result
}


#' @rdname delta_energy
#'
#' @details `ddgact_tds` calculates the entropy contribution to the change in activation energy between \code{mut} and \code{wt}
#'
#' @export
#'
ddgact_tds <- function(wt, mut, ideal = wt, pdb_site_active = NA, beta = beta_boltzmann()) {
  result <- dgact_tds(mut, ideal, pdb_site_active) - dgact_tds(wt, ideal, pdb_site_active)
  result
}



## Non-exported helper functions

#' Stress-model local-mutational-stress energy
#'
#' Calculate the energy of  `prot` in the conformation of `ideal`
#'
#' @noRd
calculate_vs <- function(prot, ideal) {
  g <- get_graph(prot)
  g_ideal <- get_graph(ideal)

  edge_in_ideal <- g$edge %in% g_ideal$edge
  g$dij_ideal <- NA
  g$dij_ideal[edge_in_ideal] <- g_ideal$dij
  g$dij_ideal[!edge_in_ideal] <- dij_edge(get_xyz(ideal), g$i[!edge_in_ideal], g$j[!edge_in_ideal])

  dij <- g$dij_ideal
  v0ij <- g$v0ij
  kij <- g$kij
  lij <- g$lij

  v_dij(dij, v0ij, kij, lij)
}
