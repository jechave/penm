#' Calculate various enm energies for prot object
#'
#' Given a protein (prot) with enm defined, it calculates energy terms
#' returns a list of internal energies (v) and entropic terms (g)
#' (note: hemholtz free energy is A = v_min + g_energy)
#'
#' @param prot A protein object, with enm graph defined
#' @param beta  boltzmann temperature
#'
#' @return A list of energy terms \code{lst(v_min,  _entropy)}
#' @export
#'
#' @examples
#'
#'@family enm builders
enm_energy <- function(prot, beta = beta_boltzmann()) {

  # internal energy terms
  v_min <- enm_v_min(prot) # energy at miniumum

  # entropic energy terms
  g_entropy <- enm_g_entropy(prot, beta)

  # return list of energy terms
  result <- lst(v_min, g_entropy)
  result
}

energy <- function(...) {
  stop("function energy renamed; call enm_energy")
}



#' Calculate minimum energy of a given prot object
#' @param prot is a prot object with a en enm$graph tibble
#' where v0ij, kij, lij and the dij for the minimum conformation are found.
#' @return a scalar: the minimum energy
enm_v_min <- function(prot) {
  graph <- prot$enm$graph
  v <- with(graph, {
    v_dij(dij, v0ij, kij, lij)
  })
  v
}

v_dij <- function(dij,v0ij,kij,lij) {
  # Calculates energy of a given conformation (dij).
  sum(v0ij + .5 * kij * (dij - lij) ^ 2)
}



#' entropic total free energy
enm_g_entropy <- function(prot, beta) {
  # Calculate T*S from the energy spectrum
  energy <- prot$enm$evalue
  sum(enm_g_entropy_mode(energy, beta))
}


#' entropic contribution of a single mode
enm_g_entropy_mode <- function(energy, beta) {
  # returns vector of entropic terms given vector of mode energies
  g_entropy_mode <- 1/(2 * beta) * log((beta * energy) / (2 * pi))
  g_entropy_mode
}
