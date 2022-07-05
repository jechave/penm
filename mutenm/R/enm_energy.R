




#' Calculate minimum energy of a given prot object
#' @param prot is a prot object, wit a component graph tibble
#' where v0ij, kij, lij and the dij for the minimum conformation are found.
#' @return a scalar: the energy at the minimum-enegy conformation
enm_v_min <- function(prot) {
  graph <- get_graph(prot)
  v <- with(graph, {
    v_dij(dij, v0ij, kij, lij)
  })
  v
}

v_dij <- function(dij, v0ij, kij, lij) {
  # Calculates energy of a given conformation (dij).
  sum(v0ij + .5 * kij * (dij - lij) ^ 2)
}



#' entropic total free energy
enm_g_entropy <- function(prot, beta) {
  # Calculate T*S from the energy spectrum
  energy <- get_evalue(prot)
  sum(enm_g_entropy_mode(energy, beta))
}


#' entropic contribution of a single mode
enm_g_entropy_mode <- function(energy, beta) {
  # returns vector of entropic terms given vector of mode energies
  g_entropy_mode <- 1 / (2 * beta) * log((beta * energy) / (2 * pi))
  g_entropy_mode
}
