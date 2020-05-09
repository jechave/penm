#' Get a single-point mutant
#'
#' Returns a mutant given a wt and a site to mutate (site_mut)
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation An integer, if 0, return \code{wt} without mutating
#' @param mut_model A string specifying mutational model ("lfenm" or "sclfenm")
#' @param mut_dl_sigma The standard deviation of a normal distribution from which edge-length perturbation is picked.
#' @param mut_sd_min An integer, only edges with \code{sdij > mut_sd_min} are mutated
#' @param seed An integer, the seed for set.seed before picking perturbations
#'
#' @return A mutated protein

#' @export
#'
#' @examples
#'
#' @family enm mutating functions
#'
get_mutant_site <- function(wt, site_mut, mutation = 0, mut_model = "lfenm", mut_dl_sigma = .3, mut_sd_min = 2, seed = 241956) {

  if (mut_model == "lfenm") {
    mut <- get_mutant_site_lfenm(wt, site_mut, mutation, mut_dl_sigma, mut_sd_min, seed)
    return(mut)
  }


  if (mut_model == "sclfenm") { # recalculate enm
    mut <- get_mutant_site_sclfenm(wt, site_mut, mutation, mut_dl_sigma, mut_sd_min, seed)
    return(mut)
  }

  stop(paste("Error in get_mutant_site,  undefined mut_model:", mut_model))

}

#' Get a single-point mutant using lfenm model
#'
#' Returns a mutant given a wt and a site to mutate (site_mut)
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation An integer, if 0, return \code{wt} without mutating
#' @param mut_dl_sigma The standard deviation of a normal distribution from which edge-length perturbation is picked.
#' @param mut_sd_min An integer, only edges with \code{sdij > mut_sd_min} are mutated
#' @param seed An integer, the seed for set.seed before picking perturbations
#'
#' @return A mutated protein

#' @export
#'
#' @examples
#'
#' @family enm mutating functions
#'
get_mutant_site_lfenm <- function(wt, site_mut, mutation, mut_dl_sigma, mut_sd_min,  seed) {

  if (mutation == 0) {
    # if mutation is 0, return wt
    return(wt)
  }

  set.seed(seed + site_mut * mutation)

  delta_lij <- get_delta_lij(wt, site_mut, mut_sd_min, mut_dl_sigma)
  dxyz <- get_dxyz(wt, delta_lij)
  mut <- wt
  mut$graph$lij = wt$graph$lij + delta_lij #TODO revise this: mut parameters are w.r.t. w0, not wt...
  mut$nodes$xyz <- wt$nodes$xyz + dxyz
  mut$graph$dij <- dij_edge(mut$nodes$xyz, mut$graph$i, mut$graph$j)
  return(mut)
}



#' Get a single-point mutant using sclfenm mode
#'
#' Returns a mutant given a wt and a site to mutate (site_mut)
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation An integer, if 0, return \code{wt} without mutating
#' @param mut_dl_sigma The standard deviation of a normal distribution from which edge-length perturbation is picked.
#' @param mut_sd_min An integer, only edges with \code{sdij > mut_sd_min} are mutated
#' @param seed An integer, the seed for set.seed before picking perturbations
#'
#' @return A mutated protein

#' @export
#'
#' @examples
#'
#' @family enm mutating functions
#'
get_mutant_site_sclfenm <- function(wt, site_mut, mutation,  mut_dl_sigma, mut_sd_min,  seed) {

  stop("ERROR: sclfenm must be fixed before this option may be run")

  if (mutation == 0) {
    # if mutation is 0, return wt
    return(wt)
  }

  set.seed(seed + site_mut * mutation)

  delta_lij <- get_delta_lij(wt, site_mut, mut_sd_min, mut_dl_sigma)
  dxyz <- get_dxyz(wt, delta_lij)
  mut <- wt
  mut$graph$lij = wt$graph$lij + delta_lij #TODO revise this: mut parameters are w.r.t. w0, not wt...
  mut$nodes$xyz <- wt$nodes$xyz + dxyz
  mut <- mutate_enm(mut)
  return(mut)

}

get_delta_lij <- function(wt, site_mut, mut_sd_min, mut_dl_sigma) {
  graph <- get_graph(wt)

  delta_lij = rep(0, nrow(get_graph(wt)))

  # pick edges to mutate

  mut_edge <- (graph$i == site_mut | graph$j == site_mut) & (graph$sdij >= mut_sd_min)
  n_mut_edge <- sum(mut_edge)
  delta_lij[mut_edge] <- rnorm(n_mut_edge, 0, mut_dl_sigma)

  delta_lij
}

get_dxyz <- function(wt, delta_lij) {
  f <- get_force(wt, delta_lij)
  cmat <- get_cmat(wt)
  nzf <- f != 0 # consider only non-zero forces, to make next step faster
  dxyz <-  crossprod(cmat[nzf, ], f[nzf]) # calculate mutant equilibrium conformation (LRA)
}


#' Get force resulting from adding delta_lij to wt
#'
#'
#' @param wt the wild-type protein
#' @param delta_lij the perturbations to the wt lij parameters
#'
#' @return A force vector of size \code{3 x nsites}
#' @export
#'
#' @examples
#' @family enm mutating functions
get_force <- function(wt, delta_lij) {

  graph <- get_graph(wt)

  stopifnot(nrow(graph) == length(delta_lij))

  i <- graph$i
  j <- graph$j
  kij <- graph$kij
  eij <- get_eij(wt)

  fij = -kij * delta_lij # Force on i in the direction from i to j.


  f <- matrix(0, nrow = 3, ncol = get_nsites(wt))
  for (k in seq(nrow(graph))) {
    ik <- i[k]
    jk <- j[k]
    f[, ik] <- f[, ik] + fij[k] * eij[k,  ]
    f[, jk] <- f[, jk] - fij[k] * eij[k,  ]
  }
  as.vector(f)
}


#' mutate enm following change in protein structure and lij parameters
#'
mutate_enm <- function(prot) {

  prot <- prot %>%
    mutate_graph() %>%
    set_enm_eij() %>%
    set_enm_kmat() %>%
    set_enm_nma()

  prot
}

#' mutate graph following change in structure
#'
mutate_graph <- function(prot) {

  # the mut graph with wt contacts but mut lij parameters...
  g1 <- get_graph(prot)

  # the "self-consistent" graph: add/rm edges according to new xyz
  g2 <- set_enm_graph(prot)$graph

  # the "frustrated" graph: keep lij for edges that haven't changed
  g2[g2$edge %in% g1$edge, "lij"] <- g1[g1$edge %in% g2$edge, "lij"]

  prot$graph <- g2
  prot
}


#' Depreciated, use mutate_enm instead
#'
enm_update <- function(...) {
  stop("enm_update has been renamed to mutate_enm")
}
