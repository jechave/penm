## ENM energies

#' Boltzmann's beta, \eqn{1 / RT}
#'
#' The default inverse temperature used by the entropic energy functions
#' ([ddg_tds()], [ddgact_tds()], [dgact_tds()]).
#'
#' @param R Boltzmann's constant per mole, in kcal/(mol K)
#' @param T absolute temperature, in Kelvin
#'
#' @returns a scalar, \eqn{1 / (RT)}, in mol/kcal
#'
#' @export
#' @family enm_energy
#'
#' @examples
#' beta_boltzmann()              # 298 K
#' beta_boltzmann(T = 310)       # body temperature
#'
beta_boltzmann <- function(R = 1.986e-3, T = 298) 1 / (R * T)


#' Calculate minimum energy of a given prot object
#'
#' Sums the elastic energy stored in every spring of the network at its current
#' conformation, \eqn{\sum_{ij} v0_{ij} + \frac{1}{2} k_{ij} (d_{ij} - l_{ij})^2}.
#' It is zero for a \code{prot} fresh from [set_enm()], where each spring's rest
#' length \code{lij} is set to its actual length \code{dij}, and positive for any
#' perturbed structure — including a wild type that is a previous generation's
#' mutant.
#'
#' @param prot is a prot object, with a component graph tibble
#' where v0ij, kij, lij and the dij for the minimum conformation are found.
#' @return a scalar: the energy at the minimum-energy conformation
#'
#' @export
#' @family enm_energy
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1, ensemble = 7)
#'
#' enm_v_min(wt)                   # 0: a network fresh from set_enm() is unstrained
#' enm_v_min(mut)                  # positive: the mutation has strained it
#'
#' # Strain accumulates along a trajectory, where each mutant becomes the wild
#' # type of the next generation. So a zero minimum energy is a property of a
#' # freshly built prot, not of "wild types".
#' p <- wt
#' for (gen in 1:3) {
#'   p <- get_mutant_site(p, site_mut = 10 * gen, mutation = 1, ensemble = 7)
#'   print(enm_v_min(p))
#' }
#'
enm_v_min <- function(prot) {
  graph <- get_graph(prot)
  v <- with(graph, {
    v_dij(dij, v0ij, kij, lij)
  })
  v
}


#' Calculate entropic total free energy of prot object
#'
#' Sums the entropic free-energy contribution of every normal mode, computed from
#' the ENM eigenvalue spectrum.
#'
#' @param prot is a prot object with known eigenvalues of enm model
#' @param beta inverse temperature, \code{1 / kT}
#'
#' @return a scalar: the entropic free-energy contribution, summed over modes
#'
#' @family enm_energy
#' @export
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' # beta has no default; beta_boltzmann() is the usual choice
#' enm_g_entropy(wt, beta_boltzmann())
#'
#' # the term depends on the whole eigenvalue spectrum, so it moves when the
#' # network does: a shorter cutoff keeps fewer contacts
#' wt2 <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'                d_max = 9)
#' enm_g_entropy(wt2, beta_boltzmann())
#'
enm_g_entropy <- function(prot, beta) {
  # Calculate T*S from the energy spectrum
  energy <- get_evalue(prot)
  sum(enm_g_entropy_mode(energy, beta))
}

## Internal

#' Calculate energy of a single spring
#'
#' @noRd
#'
v_dij <- function(dij, v0ij, kij, lij) {
  # Calculates energy of a given conformation (dij).
  sum(v0ij + .5 * kij * (dij - lij) ^ 2)
}






#' Entropic contribution of a single mode
#'
#' @noRd
#'
enm_g_entropy_mode <- function(energy, beta) {
  # returns vector of entropic terms given vector of mode energies
  g_entropy_mode <- 1 / (2 * beta) * log((beta * energy) / (2 * pi))
  g_entropy_mode
}
