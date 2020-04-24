#' Calculate various enm energies for prot object
#'
#' Given a protein (prot) with enm defined, it calculates energy terms
#' returns a list of internal energies (v) and entropic terms (g)
#' (note: hemholtz free energy is A = v_min + g_energy)
#'
#' @param prot A protein object, with enm graph defined
#' @param ideal A protein object that corresponds to the ideal "active" conformation
#' @param sd_min Integer representing a cut-off in sequence distance for energy calculations
#' @param beta  boltzmann temperature
#'
#' @return A list of energy terms \code{lst(v_min, dv_activation, g_entropy, g_entropy_activation, v_stress)}
#' @export
#'
#' @examples
#'
#'@family enm builders
enm_energy <- function(prot, ideal, sd_min = 1, beta = beta_boltzmann()) {

  # internal energy terms
  v_min <- enm_v_min(prot, sd_min) # energy at miniumum
  dv_activation <- enm_dv_activation(prot, ideal) # activation energy to shift prot$xyz to ideal$xyz for active site
  v_stress <- enm_v_stress(prot, ideal) # energy needed to shift whole protein from prot$xyz to ideal$xyz (stress model)

  # entropic energy terms
  g_entropy <- g_entropy(prot, beta)
  g_entropy_activation <- g_entropy_activation(prot, beta)

  # return list of energy terms
  lst(v_min,
      dv_activation,
      g_entropy,
      g_entropy_activation,
      v_stress)
}

#' @rdname enm_energy
#' @export
energy <- enm_energy


#' Calculate minimum energy of a given prot object
#' @param prot is a prot object with a en enm$graph tibble
#' where v0ij, kij, lij and the dij for the minimum conformation are found.
#' @return a scalar: the minimum energy
enm_v_min <- function(prot, sd_min = 1) {
  graph <- prot$enm$graph
  graph <- graph[graph$sdij >= sd_min,] # Don't include i-i+1 terms in v_min
  v_min <- with(graph, {
    v_dij(dij, v0ij, kij, lij)
  })
  v_min
}



#' Internal energy contribution to free energy of activaton
enm_dv_activation <- function(prot,ideal) {
  if (anyNA(prot$site_active)) {
    dv_activation <- NA
  } else {
    dxyz <- prot$xyz - ideal$xyz
    site_active <- prot$site_active
    dxyz <- my_as_xyz(dxyz)
    dxyz <- as.vector(dxyz[,site_active])
    kmat_active <- prot$enm$kmat_active
    dv_activation <- .5 * my_quad_form(dxyz, kmat_active, dxyz)
  }

  dv_activation
}


#' Stress-model local-mutational-stress energy
enm_v_stress <- function(prot,ideal, sd_min = 1) {
  prot_graph <- prot$enm$graph
  ideal_graph <- ideal$enm$graph

  ideal_graph <- ideal_graph  %>%
    dplyr::select(edge, dij) %>%
    rename(dij_ideal = dij)

  joint_graph <- inner_join(prot_graph, ideal_graph, by = "edge")

  #debug
  # don't include sequence neigbors in stress calculation
  joint_graph <- joint_graph %>%
    filter(sdij >= sd_min)

  kij <- joint_graph$kij
  lij <- joint_graph$lij
  v0ij <- joint_graph$v0ij
  dij <- joint_graph$dij_ideal

  v_dij(dij, v0ij, kij, lij)
}

v_dij <- function(dij,v0ij,kij,lij) {
  # Calculates energy of a given conformation (dij).
  sum(v0ij + .5 * kij * (dij - lij) ^ 2)
}



#' entropic total free energy
g_entropy <- function(prot, beta) {
  # Calculate T*S from the energy spectrum
  energy <- prot$enm$evalue
  sum(g_entropy_mode(energy, beta))
}

#' entropic contribution to dg_activation
g_entropy_activation <- function(prot, beta) {
  # Calculate entropic contribution to dg_activation
  if (anyNA(prot$enm$kmat_active)) {
    gact <- NA
  } else {
    eig <- eigen(prot$enm$kmat_active, symmetric = TRUE, only.values = TRUE)
    evalue <- eig$values
    gact <- sum(g_entropy_mode(evalue, beta))
  }
  gact

}

#' entropic contribution of a single mode
g_entropy_mode <- function(energy, beta) {
  # returns vector of entropic terms given vector of mode energies
  g_entropy <- 1/(2 * beta) * log((beta * energy) / (2 * pi))
  g_entropy
}

