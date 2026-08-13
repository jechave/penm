# Structure differences ---------------------------------------------------

#' Calculate site-dependent profiles of structure differences between two proteins
#'
#' This version works only for wt and mut with no indels
#'
#'
#' @param wt A protein object with \code{xyz} defined
#' @param mut A second protein object  with \code{xyz} defined
#'
#' @return A vector \code{(x_i)} of size \code{nsites}, where \code{x_i} is the property compared, for site i.
#'
#' @name delta_structure_by_site
#'
NULL

#' @rdname delta_structure_by_site
#'
#' @details `delta_structure_dr2i` returns the square of structural difference vector \eqn{\mathbf{C}\mathbf{f}}
#'
#' @export
#'
delta_structure_dr2i <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dxyz <- my_as_xyz(mut$nodes$xyz - wt$nodes$xyz) # use c(3, nsites) representation of xyz
  dr2i <- colSums(dxyz^2)
  dr2i
}

#' @rdname delta_structure_by_site
#' @details `delta_structure_de2i` returns the square of deformation energy vector \eqn{\mathbf{C}^{1/2}\mathbf{f}}
#'
#' @export
#'
delta_structure_de2i <- function(wt, mut, kmat_sqrt) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dr <- as.vector(get_xyz(mut) - get_xyz(wt))
  de <- kmat_sqrt %*% dr
  de <- my_as_xyz(de)
  de2i <-  colSums(de^2)
  de2i
}


#' @rdname delta_structure_by_site
#' @details `delta_structure_df2i` returns the square of force vector \eqn{\mathbf{f}}
#'
#' @export
#'
delta_structure_df2i <- function(wt, mut) {
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


#' @rdname delta_structure_by_site
#' @details `delta_structure_dvmi` returns the difference of site-dependent minimum-energy contributions
#'
#' @export
#'
delta_structure_dvmi <- function(wt, mut) {
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

#' @rdname delta_structure_by_site
#' @details `delta_structure_dvsi` returns the difference of site-dependent stress-energy contributions
#'
#' @export
#'
delta_structure_dvsi <- function(wt, mut) {
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

#' @rdname delta_structure_by_site
#' @details `delta_structure_dvsi_same_topology` returns the difference of site-dependent stress-energy contributions, assumes no change in topology
#'
#' @export
#'
delta_structure_dvsi_same_topology <- function(wt, mut) {
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


