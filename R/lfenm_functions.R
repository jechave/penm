#' Get a single-point mutant
#'
#' Returns a mutant given a wt and a site to mutate (site_mut)
#'     wt0 is the protein to use to calculate perturbations of current site
#'     (if wt0 = wt, it is the current wild-type, but it can be a fixed initial state)
#'     ideal is the ideal protein (assumed to be by default the fixed initial state)
#'     copy wt
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation An integer, if 0, return \code{wt} without mutating
#' @param wt0 A protein with respect to which mutations are to be calculated
#' @param ideal A protein with ideal protein structure
#' @param sd_min An integer, only edges with \code{sdij > sd_min} are mutated
#' @param dl_sigma The standard deviation of a normal distribution from which edge-length perturbation is picked.
#' @param model The enm model type
#' @param d_max The enm distance cut-off to define contacts
#' @param v0 The minimum energy
#' @param add_frust Whether to include frustrations or not
#' @param update_enm Whether to update enm or not
#'
#' @return
#' @export
#'
#' @examples
#'
#' @family enm mutating functions
#'
get_mutant_site <- function(wt, site_mut, mutation = 0,
                              wt0 = wt, ideal = wt0,
                              sd_min = 2, dl_sigma = .3,
                              model = "ming_wall", d_max = 10.5, v0 = 0,
                              add_frust = FALSE,
                              update_enm = FALSE) {


  if (mutation == 0) {
    # if mutation is 0, return wt
    return(wt)
  }

  mut <- wt

  # pick edges to mutate
  mut_edge <-
    (mut$enm$graph$i == site_mut |
       mut$enm$graph$j == site_mut) &
    mut$enm$graph$sdij >= sd_min

  # mutate edges
  n_mut_edge <- sum(mut_edge)
  dlij <- rnorm(n_mut_edge, 0, dl_sigma)
  #dlij <- sample(c(-dl_max,0,dl_max),size = n_mut_edge, replace = TRUE)
  mut$enm$graph$lij[mut_edge] <-  wt0$enm$graph$lij[mut_edge] + dlij

  # calculate mutant equilibrium conformation (LRA)
  f <- get_force(wt, mut)
  nzf <- f != 0 # consider only non-zero forces
  dxyz <-  crossprod(wt$enm$cmat[nzf, ], f[nzf])
  mut$xyz <- wt$xyz + dxyz


  # recalculate energies

  if (update_enm) {
    # recalculate enm
    # recalculate mutant's graph
    g1 <- mut$enm$graph

    # add/rm edges according to new xyz
    g2 <- enm_graph_xyz(mut$xyz, mut$pdb_site,
                           model = model,  d_max = d_max)

    # keep lij for edges that haven't changed
    g2[g2$edge %in% g1$edge, "lij"] <- g1[g1$edge %in% g2$edge, "lij"]

    mut$enm$graph <- g2

    # add v0ij to graph edges
    mut <- add_v0(mut, v0)

    # recalculate eij versors
    mut$enm$eij <- eij_edge(mut$xyz, mut$enm$graph$i, mut$enm$graph$j)

    # recalculate kmat
    mut$enm$kmat <-  kmat_graph(mut$enm$graph, mut$enm$eij,
                                nsites = length(mut$pdb_site),
                                add_frust = add_frust)
    # do nma and update nma info in mut
    nma <- enm_nma(mut$enm)
    mut$enm$mode = nma$mode
    mut$enm$evalue = nma$evalue
    mut$enm$umat = nma$umat
    mut$enm$cmat = nma$cmat

    # if site_active is defined, recalculate cmat_active and kmat_active
    if (anyNA(mut$ind_active)) {
      mut$enm$cmat_active <- NA
      mut$enm$kmat_active <- NA
    } else {
      mut$enm$cmat_active <- mut$enm$cmat[mut$ind_active, mut$ind_active]
      mut$enm$kmat_active <- solve(mut$enm$cmat_active)
    }

  } else {
  # recalculate only equilibrium distances
    mut$enm$graph$dij <- dij_edge(mut$xyz, mut$enm$graph$i, mut$enm$graph$j)
  }


  # calculate mutant energy
  mut$energy <- energy(mut, ideal = ideal, sd_min = 1)

  mut
}



#' Get force that mutates wt into mut
#'
#' Given \code{mut} obtained by perturbing site edges of \code{wt}, obtain force \code{fij = kij * dlij}
#'
#' @param wt The wild-type protein
#' @param mut The mutant, which contains the mutated \code{graph}
#'
#' @return A force vector of size \code{3 x nsites}
#' @export
#'
#' @examples
#' @family enm mutating functions
get_force <- function(wt, mut) {
  n_edges <- nrow(wt$enm$graph)
  nsites <- wt$nsites
  i <- wt$enm$graph$i
  j <- wt$enm$graph$j
  dij <- wt$enm$graph$dij # equilibrium of wt
  kij <- wt$enm$graph$kij
  lij <- mut$enm$graph$lij # mutant's spring lengths
  fij = kij*(dij - lij) # Force on i in the direction from i to j.

  eij <- wt$enm$eij
  f <- matrix(0, nrow = 3, ncol = nsites)
  for (k in seq(n_edges)) {
    ik <- i[k]
    jk <- j[k]
    f[, ik] <- f[, ik] + fij[k] * eij[k,  ]
    f[, jk] <- f[, jk] - fij[k] * eij[k,  ]
  }
  as.vector(f)
}

#' Update ENM following graph change
#'
#' @param prot  A protein object
#' @param add_frust A logical indicating whether to add frustration
#'
#' @return A protein object with updated \code{enm} component
#' @export
#'
#' @examples
#'
#' @family enm mutating functions
#'
enm_update <- function(prot, add_frust = FALSE) {
  # it re-calculates what needs to be recalculated due to change in graph
  stopifnot(!is.null(prot$ind_active)) # stop if no active-diste info
  enm <- prot$enm
  enm$kmat <- kmat_graph(prot$enm$graph, prot$enm$eij, prot$nsites, add_frust)
  nma <- enm_nma(enm) #returns mode, evalue, cmat, umat
  enm$mode <- nma$mode
  enm$evalue <- nma$evalue
  enm$cmat <- nma$cmat
  enm$umat <- nma$umat

  if (anyNA(prot$ind_active)) { #if ind_active is NA, cmat_active and kmat_active are undefined
    enm$cmat_active <- NA
    enm$kmat_active <- NA
  } else {
    enm$cmat_active <- enm$cmat[prot$ind_active, prot$ind_active]
    enm$kmat_active <- solve(enm$cmat_active)
  }

  prot$enm <- enm
  prot
}

#' @rdname enm_update
#' @export
update_enm <- enm_update

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
  stopifnot(prot1$pdb_site == prot2$pdb_site) # this version of dr2_site is to compare proteins with no indels
  stopifnot(prot1$site == prot2$site) # this version of dr2_site is to compare proteins with no indels
  dxyz <- my_as_xyz(prot2$xyz - prot1$xyz) # use c(3, nsites) representation of xyz
  dr2i <- colSums(dxyz^2)
  dr2i
}

#' @rdname dr2_site
de2_site <- function(prot1, prot2) {
  stopifnot(prot1$pdb_site == prot2$pdb_site) # this version of dr2_site is to compare proteins with no indels
  stopifnot(prot1$site == prot2$site) # this version of dr2_site is to compare proteins with no indels

  evalue <- prot1$enm$evalue
  umat <- prot1$enm$umat
  kmat_sqrt <- umat %*% (sqrt(evalue) * t(umat))
  dxyz <- my_as_xyz(prot2$xyz - prot1$xyz) # use c(3, nsites) representation of xyz
  dr <- as.vector(dxyz)
  de <- kmat_sqrt %*% dr
  de <- my_as_xyz(de)
  de2i <-  colSums(de^2)
  de2i
}

#' @rdname dr2_site
df2_site <- function(prot1, prot2) {
  stopifnot(prot1$pdb_site == prot2$pdb_site) # this version of dr2_site is to compare proteins with no indels
  stopifnot(prot1$site == prot2$site) # this version of dr2_site is to compare proteins with no indels

  kmat <- prot1$enm$kmat

  dxyz <- my_as_xyz(prot2$xyz - prot1$xyz) # use c(3, nsites) representation of xyz
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
  stopifnot(prot1$pdb_site == prot2$pdb_site) # this version of dr2_site is to compare proteins with no indels
  dxyz <- as.vector(prot2$xyz - prot1$xyz) # use c(3, nsites) representation of xyz
  drn <- as.vector(crossprod(prot1$enm$umat, dxyz))
  stopifnot(length(prot1$enm$mode) == length(drn))
  dr2n <- drn^2
  as.vector(dr2n)
}

#' @rdname dr2_nm
de2_nm <- function(prot1, prot2) {
  stopifnot(prot1$pdb_site == prot2$pdb_site) # this version of dr2_site is to compare proteins with no indels
  stopifnot(prot1$site == prot2$site) # this version of dr2_site is to compare proteins with no indels

  evalue <- prot1$enm$evalue
  umat <- prot1$enm$umat
  kmat_sqrt <- umat %*% (sqrt(evalue) * t(umat))
  dr <- as.vector(prot2$xyz - prot1$xyz) # use c(3, nsites) representation of xyz
  de <- kmat_sqrt %*% dr
  den <- t(umat) %*% de
  de2n <-  den^2
  as.vector(de2n)
}



#' @rdname dr2_nm
df2_nm <- function(prot1, prot2) {
  stopifnot(prot1$pdb_site == prot2$pdb_site) # this version of dr2_site is to compare proteins with no indels

  kmat <- prot1$enm$kmat
  umat <- prot1$enm$umat

  dr <- as.vector(prot2$xyz - prot1$xyz)
  df <- kmat %*% dr
  dfn <- t(umat) %*% df
  df2n <-  dfn^2
  as.vector(df2n)
}






