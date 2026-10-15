## ENM energies

#' Boltzmann's beta, \eqn{1 / RT}{1 / RT}
#'
#' The default inverse temperature used by the entropic energy functions
#' ([ddg_tds()], [ddgact_tds()], [dgact_tds()]).
#'
#' penm works in kcal/mol: every model's spring constants are in
#' kcal/(mol Å²) (see [set_enm()]), so its energies are in kcal/mol and the
#' default `R`, in kcal/(mol K), matches them.
#'
#' @param R the gas constant, in kcal/(mol K); the default is the CODATA 2018
#'   value, 8.314462618 J/(mol K), divided by 4184 J/kcal
#' @param T absolute temperature, in Kelvin
#'
#' @returns a scalar, \eqn{1 / (RT)}{1 / (RT)}, in mol/kcal
#'
#' @export
#' @family enm_energy
#'
#' @examples
#' beta_boltzmann()              # 298 K
#' beta_boltzmann(T = 310)       # body temperature
#'
beta_boltzmann <- function(R = 1.987204e-3, T = 298) 1 / (R * T)


#' Calculate minimum energy of a given prot object
#'
#' Sums the elastic energy stored in every spring of the network at its current
#' conformation, \eqn{\sum_{ij} v0_{ij} + \frac{1}{2} k_{ij} (d_{ij} - l_{ij})^2}{sum_ij [v0_ij + (1/2) k_ij (d_ij - l_ij)^2]}.
#' It is zero for a \code{prot} fresh from [set_enm()], where each spring's rest
#' length \code{lij} is set to its actual length \code{dij}, and positive for any
#' perturbed structure — including a wild type that is a previous generation's
#' mutant.
#'
#' @param prot is a prot object, with a component graph tibble
#' where v0ij, kij, lij and the dij for the minimum conformation are found.
#' @return a scalar: the energy at the minimum-energy conformation, in kcal/mol
#'
#' @export
#' @family enm_energy
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, ensemble = 7)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1)
#'
#' enm_v_min(wt)                   # 0: a network fresh from set_enm() is unstrained
#' enm_v_min(mut)                  # positive: the mutation has strained it
#'
#' # Strain accumulates along a trajectory, where each mutant becomes the wild
#' # type of the next generation. So a zero minimum energy is a property of a
#' # freshly built prot, not of "wild types".
#' p <- wt
#' for (gen in 1:3) {
#'   p <- get_mutant_site(p, site_mut = 10 * gen, mutation = 1)
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
#' the ENM eigenvalue spectrum, \eqn{\frac{1}{2\beta} \sum_n \ln \frac{\beta\lambda_n}{2\pi}}{sum_n (1/(2 beta)) ln(beta lambda_n / (2 pi))}.
#' Its value depends on the unit of length (Å); differences between two networks
#' with the same number of modes ([ddg_tds()]) do not.
#'
#' @param prot is a prot object with known eigenvalues of enm model
#' @param beta inverse temperature \eqn{1 / RT}{1 / RT}, in mol/kcal (see [beta_boltzmann()])
#'
#' @return a scalar: the entropic free-energy contribution, summed over modes, in
#'   kcal/mol
#'
#' @family enm_energy
#' @export
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' enm_g_entropy(wt)                           # at 298 K
#' enm_g_entropy(wt, beta_boltzmann(T = 310))  # at 310 K
#'
#' # the term depends on the whole eigenvalue spectrum, so it moves when the
#' # network does: a shorter cutoff keeps fewer contacts
#' wt2 <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'                d_max = 9)
#' enm_g_entropy(wt2, beta_boltzmann())
#'
enm_g_entropy <- function(prot, beta = beta_boltzmann()) {
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
