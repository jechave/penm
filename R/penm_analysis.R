calculate_dvm <- function(prot1, prot2)
  enm_v_min(prot2) - enm_v_min(prot1)

delta_g_entropy <- function(prot1, prot2, beta)
  enm_g_entropy(prot2, beta) - enm_g_entropy(prot1, beta)


#' Stress-model difference of local-mutational-stress energy
calculate_dvs <- function(prot1, prot2, ideal = prot1)
  enm_v_stress(prot2, ideal) - enm_v_stress(prot1, ideal)


#' Stress-model local-mutational-stress energy
enm_v_stress <- function(prot, ideal) {
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
#' This version works only for prot1 and prot2 with no indels
#'
#' @param prot1 A protein object with \code{xyz} defined
#' @param prot2 A second protein object  with \code{xyz} defined
#'
#' @return A tibble with columns \code{pdb_site, site, dr2}
#' @export
#'
#' @examples
calculate_dr2i <- function(prot1, prot2) {
  stopifnot(prot1$node$pdb_site == prot2$node$pdb_site) # no indels
  stopifnot(prot1$node$site == prot2$node$site) # no indels
  dxyz <- my_as_xyz(prot2$nodes$xyz - prot1$nodes$xyz) # use c(3, nsites) representation of xyz
  dr2i <- colSums(dxyz^2)
  dr2i
}

#' @rdname calculate_dr2i
calculate_de2i <- function(prot1, prot2, kmat_sqrt) {
  stopifnot(prot1$node$pdb_site == prot2$node$pdb_site) # no indels
  stopifnot(prot1$node$site == prot2$node$site) # no indels
  dr <- as.vector(get_xyz(prot2) - get_xyz(prot1))
  de <- kmat_sqrt %*% dr
  de <- my_as_xyz(de)
  de2i <-  colSums(de^2)
  de2i
}



#' @rdname calculate_dr2i
calculate_df2i <- function(prot1, prot2) {
  stopifnot(prot1$node$pdb_site == prot2$node$pdb_site) # no indels
  stopifnot(prot1$node$site == prot2$node$site) # no indels

  kmat <- prot1$kmat

  dxyz <- my_as_xyz(prot2$nodes$xyz - prot1$nodes$xyz) # use c(3, nsites) representation of xyz
  dr <- as.vector(dxyz)
  df <- kmat %*% dr
  df <- my_as_xyz(df)
  df2i <-  colSums(df^2)
  df2i
}

dvm_site <- function(prot1, prot2) {
  stopifnot(get_nsites(prot1) == get_nsites(prot2)) # #warning, #check: I'm assuming no indels

  g1 <- get_graph(prot1) %>%
    mutate(vmij = v0ij + 1/2 * kij * (dij - lij)^2)
  g2 <- get_graph(prot2) %>%
    mutate(vmij = v0ij + 1/2 * kij * (dij - lij)^2)

  g <- inner_join(g1, g2, by = c("edge", "i", "j")) %>%
    mutate(dvmij = vmij.y - vmij.x)


  nsites <- get_nsites(prot1)
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
#' This version works only for prot1 and prot2 with no indels
#'
#' @param prot1 A protein object with \code{xyz} and \code{enm} defined
#' @param prot2 A second protein object  with \code{xyz} defined
#'
#' @return A vector of square differences along normal modes \code{dr2n}
#' @export
#'
#' @examples
calculate_dr2n <- function(prot1, prot2) {
  stopifnot(prot1$node$pdb_site == prot2$node$pdb_site) # no indels
  dr <- as.vector(get_xyz(prot2) - get_xyz(prot1))
  drn <- as.vector(crossprod(get_umat(prot1), dr))
  stopifnot(length(prot1$nma$mode) == length(drn))
  dr2n <- drn^2
  as.vector(dr2n)
}


#' @rdname calculate_dr2n
calculate_de2n <- function(prot1, prot2) {
  get_evalue(prot1) * calculate_dr2n(prot1, prot2)
}

#' @rdname calculate_dr2n
calculate_df2n <- function(prot1, prot2) {
  get_evalue(prot1)^2 * calculate_dr2n(prot1, prot2)
}
