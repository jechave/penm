#' add active-site indexes to prot object
#'
add_active_site_indexes <- function(prot, pdb_site_active) {
  result <- prot
  result$pdb_site_active <- pdb_site_active # add pdb_site_active
  result$site_active <- prot$site[prot$pdb_site %in% pdb_site_active] # add site_active
  result$ind_active <- site_to_ind(result$site_active) # add ind_active
  result
}

#' Add cmat_activation and kmat_activation to prot object
#'
#' @param prot A protein object that must contain \code{xyz} and \code{pdb_site} elements, and active_site indexes \code{ind_active} (which may be NA)
#'
#' @return A protein object with added cmat_active and kmat_active
#'
#' @export
#'
#' @family enm builders
#'
#' @examples
add_enm_activation <- function(prot)  {
  result <- prot

    result$enm$cmat_active <- result$enm$cmat[result$ind_active, result$ind_active]
    result$enm$kmat_active <- solve(result$enm$cmat_active)

  result
}

mutate_enm_activation <- add_enm_activation




#' Stress-model local-mutational-stress energy
enm_v_stress <- function(prot, ideal) {
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
enm_dv_activation <- function(prot, ideal) {
  if (anyNA(prot$site_active)) {
    dv_activation <- NA
  } else {
    dxyz <- prot$xyz - ideal$xyz
    site_active <- prot$site_active
    dxyz <- my_as_xyz(dxyz)
    dxyz <- as.vector(dxyz[, site_active])
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


distance_to_active <- function(xyz, site_active) {
  stopifnot(length(xyz) %% 3 == 0)
  nsites <- length(xyz) / 3

  if (anyNA(site_active)) {
    distance <-  rep(NA, nsites)
  } else {
    site <- seq(nsites)
    is_active_site <- site %in% site_active
    xyz <- my_as_xyz(xyz)

    nsite_active <-  sum(is_active_site)
    distance <- rep(NA, nsites)
    for (j in seq(nsites)) {
      d_active_to_j <-  xyz[, j] - xyz[, is_active_site]
      dim(d_active_to_j) <-  c(3, nsite_active)
      d_active_to_j_norm <-  sqrt(colSums(d_active_to_j^2))
      distance[j] <-  min(d_active_to_j_norm)
    }
  }

  distance
}


delta_v_stress <- function(prot1, prot2)
  enm_v_stress(prot2) - enm_v_stress(prot1)

delta_v_activation <- function(prot1, prot2)
  enm_dv_activation(prot2) - enm_dv_activation(prot1)

delta_g_entropy_activation <- function(prot1, prot2)
  enm_g_entropy_activation(prot2) - enm_g_entropy_activation(prot1)
