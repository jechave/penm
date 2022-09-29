## Energy diferences
#' Calculate energy differences between a mutant and wild type
#'
#' @param wt A protein object with \code{xyz} defined
#' @param mut A second protein object  with \code{xyz} defined
#'
#' @return A (scalar) energy difference between mutant and wild type.
#'
#' @name delta_energy
#'
NULL

#' @rdname delta_energy
#'
#' @details `calculate_dvm` calculates the minimum-energy difference between \code{mut} and \code{wt}
#'
#' @export
#'
calculate_dvm <- function(wt, mut)
  enm_v_min(mut) - enm_v_min(wt)

#' @rdname delta_energy
#'
#' @details `calculate_dg_enetropy` calculates the entropic free energ difference between \code{mut} and \code{wt}
#'
#' @export
#'
calculate_dg_entropy <- function(wt, mut, beta)
  enm_g_entropy(mut, beta) - enm_g_entropy(wt, beta)


#' @rdname delta_energy
#'
#' @details `calculate_dvs` calculates the ideal-conformation stress-energy difference between \code{mut} and \code{wt}
#'
#' @export
#'
calculate_dvs <- function(wt, mut, ideal = wt)
  calculate_vs(mut, ideal) - calculate_vs(wt, ideal)


#' Stress-model local-mutational-stress energy
#'
#' Calculate de energy of a configuration "ideal"
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
