delta_v_min <- function(prot1, prot2)
  enm_v_min(prot2) - enm_v_min(prot1)

delta_g_entropy <- function(prot1, prot2)
  enm_g_entropy(prot2) - enm_g_entropy(prot1)




#' Compare the structures of two proteins in site representation
#'
#' Calculate de difference between the structures of two proteins, return dr2_site.
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
dr2_site <- function(prot1, prot2) {
  stopifnot(prot1$node$pdb_site == prot2$node$pdb_site) # this version of dr2_site is to compare proteins with no indels
  stopifnot(prot1$node$site == prot2$node$site) # this version of dr2_site is to compare proteins with no indels
  dxyz <- my_as_xyz(prot2$nodes$xyz - prot1$nodes$xyz) # use c(3, nsites) representation of xyz
  dr2i <- colSums(dxyz^2)
  dr2i
}

#' @rdname dr2_site
de2_site <- function(prot1, prot2) {
  stopifnot(prot1$node$pdb_site == prot2$node$pdb_site) # this version of dr2_site is to compare proteins with no indels
  stopifnot(prot1$node$site == prot2$node$site) # this version of dr2_site is to compare proteins with no indels

  evalue <- prot1$nma$evalue
  umat <- prot1$nma$umat
  kmat_sqrt <- umat %*% (sqrt(evalue) * t(umat))
  dxyz <- my_as_xyz(prot2$nodes$xyz - prot1$nodes$xyz) # use c(3, nsites) representation of xyz
  dr <- as.vector(dxyz)
  de <- kmat_sqrt %*% dr
  de <- my_as_xyz(de)
  de2i <-  colSums(de^2)
  de2i
}

#' @rdname dr2_site
df2_site <- function(prot1, prot2) {
  stopifnot(prot1$node$pdb_site == prot2$node$pdb_site) # this version of dr2_site is to compare proteins with no indels
  stopifnot(prot1$node$site == prot2$node$site) # this version of dr2_site is to compare proteins with no indels

  kmat <- prot1$kmat

  dxyz <- my_as_xyz(prot2$nodes$xyz - prot1$nodes$xyz) # use c(3, nsites) representation of xyz
  dr <- as.vector(dxyz)
  df <- kmat %*% dr
  df <- my_as_xyz(df)
  df2i <-  colSums(df^2)
  df2i
}



#' Compare the structures of two proteins in nm representation
#'
#' Calculate de difference between the structures of two proteins, return dr2_nm.
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
dr2_nm <- function(prot1, prot2) {
  stopifnot(prot1$node$pdb_site == prot2$node$pdb_site) # this version of dr2_site is to compare proteins with no indels
  dxyz <- as.vector(prot2$nodes$xyz - prot1$nodes$xyz) # use c(3, nsites) representation of xyz
  drn <- as.vector(crossprod(prot1$nma$umat, dxyz))
  stopifnot(length(prot1$nma$mode) == length(drn))
  dr2n <- drn^2
  as.vector(dr2n)
}

#' @rdname dr2_nm
de2_nm <- function(prot1, prot2) {
  stopifnot(prot1$node$pdb_site == prot2$node$pdb_site) # this version of dr2_site is to compare proteins with no indels
  stopifnot(prot1$node$site == prot2$node$site) # this version of dr2_site is to compare proteins with no indels

  evalue <- prot1$nma$evalue
  umat <- prot1$nma$umat
  kmat_sqrt <- umat %*% (sqrt(evalue) * t(umat))
  dr <- as.vector(prot2$nodes$xyz - prot1$nodes$xyz) # use c(3, nsites) representation of xyz
  de <- kmat_sqrt %*% dr
  den <- t(umat) %*% de
  de2n <-  den^2
  as.vector(de2n)
}



#' @rdname dr2_nm
df2_nm <- function(prot1, prot2) {
  stopifnot(prot1$node$pdb_site == prot2$node$pdb_site) # this version of dr2_site is to compare proteins with no indels

  kmat <- prot1$kmat
  umat <- prot1$nma$umat

  dr <- as.vector(prot2$nodes$xyz - prot1$nodes$xyz)
  df <- kmat %*% dr
  dfn <- t(umat) %*% df
  df2n <-  dfn^2
  as.vector(df2n)
}
