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



# Structure differences ---------------------------------------------------
#' Compare two protein structures site by site (site-dependent profiles)
#'
#' This version works only for wt and mut with no indels
#'
#'
#' @param wt A protein object with \code{xyz} defined
#' @param mut A second protein object  with \code{xyz} defined
#'
#' @return A vector \code{(x_i)} of size \code{nsites}, where \code{x_i} is the property compared, for site i.
#'
#' @name site_profile
#'
NULL

#' @rdname site_profile
#'
#' @details `calculate_dr2i` returns the square of structural difference vector \eqn{\mathbf{C}\mathbf{f}}
#'
#' @export
#'
calculate_dr2i <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dxyz <- my_as_xyz(mut$nodes$xyz - wt$nodes$xyz) # use c(3, nsites) representation of xyz
  dr2i <- colSums(dxyz^2)
  dr2i
}

#' @rdname site_profile
#' @details `calculate_de2i` returns the square of deformation energy vector \eqn{\mathbf{C}^{1/2}\mathbf{f}}
#'
#' @export
#'
calculate_de2i <- function(wt, mut, kmat_sqrt) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dr <- as.vector(get_xyz(mut) - get_xyz(wt))
  de <- kmat_sqrt %*% dr
  de <- my_as_xyz(de)
  de2i <-  colSums(de^2)
  de2i
}


#' @rdname site_profile
#' @details `calculate_df2i` returns the square of force vector \eqn{\mathbf{f}}
#'
#' @export
#'
calculate_df2i <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels

  kmat <- wt$kmat

  dxyz <- my_as_xyz(mut$nodes$xyz - wt$nodes$xyz) # use c(3, nsites) representation of xyz
  dr <- as.vector(dxyz)
  df <- kmat %*% dr
  df <- my_as_xyz(df)
  df2i <-  colSums(df^2)
  df2i
}


#' @rdname site_profile
#' @details `calculate_dvmi` returns the difference of site-dependent minimum-energy contributions
#'
#' @export
#'
calculate_dvmi <- function(wt, mut) {
  stopifnot(get_nsites(wt) == get_nsites(mut)) # #warning, #check: I'm assuming no indels

  g1 <- get_graph(wt) %>%
    mutate(vmij = v0ij + 1/2 * kij * (dij - lij)^2)
  g2 <- get_graph(mut) %>%
    mutate(vmij = v0ij + 1/2 * kij * (dij - lij)^2)

  g <- inner_join(g1, g2, by = c("edge", "i", "j")) %>%
    mutate(dvmij = vmij.y - vmij.x)


  nsites <- get_nsites(wt)
  dvmi <- rep(0, nsites)

  for (site in seq(nsites)) {
    dvmi[site] <-  sum(g$dvmij[g$i == site | g$j == site])
  }

  dvmi
}

#' @rdname site_profile
#' @details `calculate_dvsi` returns the difference of site-dependent stress-energy contributions
#'
#' @export
#'
calculate_dvsi <- function(wt, mut) {
  stopifnot(get_nsites(wt) == get_nsites(mut)) # #warning, #check: I'm assuming no indels

  gwt <- get_graph(wt)
  gmut <- get_graph(mut)

  stopifnot(all(gmut$edge == gwt$edge)) # for "lfenm": this works if the network didn't change its topology

  g <- inner_join(gwt, gmut, by = c("edge", "i", "j"), suffix = c(".wt", ".mut"))  %>%
    mutate(dvsij = 1/2 * kij.mut * (dij.wt - lij.mut)^2 - 1/2 * kij.wt * (dij.wt - lij.wt)^2)

  nsites <- get_nsites(wt)
  dvsi <- rep(0, nsites)

  for (site in seq(nsites)) {
    dvsi[site] <-  sum(g$dvsij[g$i == site | g$j == site])
  }
  dvsi
}

#' @rdname site_profile
#' @details `calculate_dvsi_same_topology` returns the difference of site-dependent stress-energy contributions, assumes no change in topology
#'
#' @export
#'
calculate_dvsi_same_topology <- function(wt, mut) {
  gwt <- get_graph(wt)
  gmut <- get_graph(mut)

  stopifnot(all(gmut$edge == gwt$edge)) # for "lfenm": this works if the network didn't change its topology

  gmut$vsij = 1/2 * gmut$kij * (gwt$dij - gmut$lij)^2
  gwt$vsij = 1/2 * gwt$kij * (gwt$dij - gwt$lij)^2
  dvsij = gmut$vsij - gwt$vsij

  dvsij_non_zero <- !near(dvsij, 0)

  dvsij <- dvsij[dvsij_non_zero]
  i <- gwt$i[dvsij_non_zero]
  j <- gwt$j[dvsij_non_zero]
  sites_non_zero <- unique(c(i,j))

  dvsi <- rep(0, get_nsites(wt))

  for (e in seq_along(dvsij))  {
    dvsi[i[e]] = dvsi[i[e]] + dvsij[e]
    dvsi[j[e]] = dvsi[j[e]] + dvsij[e]
  }

  dvsi
}



#' Compare two protein structures in nm representation
#'
#' (This version works only for wt and mut with no indels)
#'
#' @param wt A protein object with \code{xyz} and \code{enm} defined
#' @param mut A second protein object  with \code{xyz} defined
#' @return A vector with contributions of each normal mode to the given property
#'
#' @name mode_profile
#'
NULL


#' @rdname mode_profile
#'
#' @details `calculate_dr2n` calculates de square of the mode-contributions to \eqn{\delta \mathbf{r} = \mathbf{C}\mathbf{f}}
#'
#'
#' @export
#'
calculate_dr2n <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dr <- as.vector(get_xyz(mut) - get_xyz(wt))
  drn <- as.vector(crossprod(get_umat(wt), dr))
  stopifnot(length(wt$nma$mode) == length(drn))
  dr2n <- drn^2
  as.vector(dr2n)
}

#' @rdname mode_profile
#'
#' @details `calculate_de2n` calculates de square of the mode-contributions to \eqn{\delta \mathbf{e} = \mathbf{C}^{1/2}\mathbf{f}}
#'
#'
#' @export
#'
calculate_de2n <- function(wt, mut) {
  get_evalue(wt) * calculate_dr2n(wt, mut)
}

#' @rdname mode_profile
#'
#' @details `calculate_df2n` calculates de square of the mode-contributions to the force vecgtor \eqn{\mathbf{f}}
#'
#'
#' @export
#'
calculate_df2n <- function(wt, mut) {
  get_evalue(wt)^2 * calculate_dr2n(wt, mut)
}


# Compare protein ensembles ----------------------------------------------------

#' Compare two protein nsembles site by site (site-dependent profiles)
#'
#' Given two proteins, compare their conformational ensembles (fluctuation patterns),
#' the proteins'  \code{cmat}, and normal modes (Principal components) are assumed known.
#'
#' (This version works only for wt and mut with no indels)
#'
#' @param wt A protein object with \code{enm} defined
#' @param mut A second protein object  with \code{enm} defined
#'
#' @return A vector \code{(x_i)} of size \code{nsites}, where \code{x_i} is the property compared, for site i.
#'
#' @name site_ensemble_profile
#'
NULL

#' @rdname site_ensemble_profile
#'
#' @details `calculate_dmsfi` returns change profile of changes of mean-square fluctuations \eqn{\delta \sigma_i^2}
#'
#' @export
#'
calculate_dmsfi <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels


}
