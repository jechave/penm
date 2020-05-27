calculate_dvm <- function(wt, mut)
  enm_v_min(mut) - enm_v_min(wt)

calculate_dg_entropy <- function(wt, mut, beta)
  enm_g_entropy(mut, beta) - enm_g_entropy(wt, beta)


#' Stress-model difference of local-mutational-stress energy
calculate_dvs <- function(wt, mut, ideal = wt)
  calculate_vs(mut, ideal) - calculate_vs(wt, ideal)


#' Stress-model local-mutational-stress energy
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



#' Compare the structures of two proteins in site representation
#'
#' Calculate de difference between the structures of two proteins, return dr2i
#'
#' This version works only for wt and mut with no indels
#'
#' @param wt A protein object with \code{xyz} defined
#' @param mut A second protein object  with \code{xyz} defined
#'
#' @return A tibble with columns \code{pdb_site, site, dr2}
#' @export
#'
#' @examples
calculate_dr2i <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dxyz <- my_as_xyz(mut$nodes$xyz - wt$nodes$xyz) # use c(3, nsites) representation of xyz
  dr2i <- colSums(dxyz^2)
  dr2i
}

#' @rdname calculate_dr2i
calculate_de2i <- function(wt, mut, kmat_sqrt) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  stopifnot(wt$node$site == mut$node$site) # no indels
  dr <- as.vector(get_xyz(mut) - get_xyz(wt))
  de <- kmat_sqrt %*% dr
  de <- my_as_xyz(de)
  de2i <-  colSums(de^2)
  de2i
}



#' @rdname calculate_dr2i
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



#' Compare the structures of two proteins in nm representation
#'
#' Calculate de difference between the structures of two proteins, return dr2n.
#'
#' This version works only for wt and mut with no indels
#'
#' @param wt A protein object with \code{xyz} and \code{enm} defined
#' @param mut A second protein object  with \code{xyz} defined
#'
#' @return A vector of square differences along normal modes \code{dr2n}
#' @export
#'
#' @examples
calculate_dr2n <- function(wt, mut) {
  stopifnot(wt$node$pdb_site == mut$node$pdb_site) # no indels
  dr <- as.vector(get_xyz(mut) - get_xyz(wt))
  drn <- as.vector(crossprod(get_umat(wt), dr))
  stopifnot(length(wt$nma$mode) == length(drn))
  dr2n <- drn^2
  as.vector(dr2n)
}


#' @rdname calculate_dr2n
calculate_de2n <- function(wt, mut) {
  get_evalue(wt) * calculate_dr2n(wt, mut)
}

#' @rdname calculate_dr2n
calculate_df2n <- function(wt, mut) {
  get_evalue(wt)^2 * calculate_dr2n(wt, mut)
}
