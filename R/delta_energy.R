## Energy diferences
#' Calculate energy differences between a mutant and wild type
#'
#' @param wt A protein object
#' @param mut A second protein
#'
#' @return A (scalar) energy difference between mutant and wild type.
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
