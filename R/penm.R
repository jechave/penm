#' Get a single-point mutant
#'
#' Returns a mutant given a wt and a site to mutate (site_mut)
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation An integer, if 0, return \code{wt} without mutating
#' @param mut_sd_min An integer, only edges with \code{sdij > mut_sd_min} are mutated
#' @param dl_sigma The standard deviation of a normal distribution from which edge-length perturbation is picked.
#' @param update_enm A logical, whether to recalculate kmat and nma or not.
#' @param seed An integer, the seed for set.seed before picking perturbations
#'
#' @return A mutated protein

#' @export
#'
#' @examples
#'
#' @family enm mutating functions
#'
get_mutant_site <- function(wt, site_mut, mutation = 0,
                            mut_sd_min = 2, dl_sigma = .3, update_enm = FALSE,
                            seed = 241956) {

  if (mutation == 0) {
    # if mutation is 0, return wt
    return(wt)
  }

  set.seed(seed + site_mut * mutation)


  delta_lij <- get_delta_lij(wt, site_mut, mut_sd_min, dl_sigma)
  f <- get_force(wt, delta_lij)
  cmat <- get_cmat(wt)
  nzf <- f != 0 # consider only non-zero forces, to make next step faster
  dxyz <-  crossprod(cmat[nzf, ], f[nzf]) # calculate mutant equilibrium conformation (LRA)



  mut <- wt
  mut$graph$lij = wt$graph$lij + delta_lij #TODO revise this: mut parameters are w.r.t. w0, not wt...

  mut$nodes$xyz <- wt$nodes$xyz + dxyz

  # recalculate enm

  if (update_enm) { # recalculate enm
    # recalculate (graph, kmat, etc. etc.)
    mut <- mutate_enm(mut)
  } else {
    # recalculate only equilibrium distances
    mut$graph$dij <- dij_edge(mut$nodes$xyz, mut$graph$i, mut$graph$j)
  }

  mut
}

get_delta_lij <- function(wt, site_mut, mut_sd_min, dl_sigma) {
  graph <- get_graph(wt)

  delta_lij = rep(0, nrow(get_graph(wt)))

  # pick edges to mutate

  mut_edge <- (graph$i == site_mut | graph$j == site_mut) & (graph$sdij >= mut_sd_min)
  n_mut_edge <- sum(mut_edge)
  delta_lij[mut_edge] <- rnorm(n_mut_edge, 0, dl_sigma)

  delta_lij
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
