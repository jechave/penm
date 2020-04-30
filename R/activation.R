# All functions related to activation, either of whole protein (v_stress) or of active site

#' Stress-model local-mutational-stress energy
enm_v_stress <- function(prot,ideal) {
  prot_graph <- prot$enm$graph
  ideal_graph <- ideal$enm$graph

  ideal_graph <- ideal_graph  %>%
    dplyr::select(edge, dij) %>%
    rename(dij_ideal = dij)

  joint_graph <- inner_join(prot_graph, ideal_graph, by = "edge")

  kij <- joint_graph$kij
  lij <- joint_graph$lij
  v0ij <- joint_graph$v0ij
  dij <- joint_graph$dij_ideal

  v_dij(dij, v0ij, kij, lij)
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

#' entropic contribution to dg_activation
enm_g_entropy_activation <- function(prot, beta) {
  # Calculate entropic contribution to dg_activation
  if (anyNA(prot$enm$kmat_active)) {
    gact <- NA
  } else {
    eig <- eigen(prot$enm$kmat_active, symmetric = TRUE, only.values = TRUE)
    evalue <- eig$values
    gact <- sum(enm_g_entropy_mode(evalue, beta))
  }
  gact

}
