#' Get a single-point mutant
#'
#' Returns a mutant given a wt and a site to mutate (site_mut)
#'     wt0 is the protein to use to calculate perturbations of current site
#'     (if wt0 = wt, it is the current wild-type, but it can be a fixed initial state)
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation An integer, if 0, return \code{wt} without mutating
#' @param wt0 A protein with respect to which mutations are to be calculated
#' @param mut_sd_min An integer, only edges with \code{sdij > mut_sd_min} are mutated
#' @param dl_sigma The standard deviation of a normal distribution from which edge-length perturbation is picked.
#' @param seed An integer, the seed for set.seed before picking perturbations
#'
#' @return A mutated protein

#' @export
#'
#' @examples
#'
#' @family enm mutating functions
get_mutant_site <- function(wt, site_mut, mutation = 0,
                            seed = 241956 + site_mut * mutation,
                            wt0 = wt,
                            mut_sd_min = 2, dl_sigma = .3, update_enm = FALSE) {

  if (mutation == 0) {
    # if mutation is 0, return wt
    return(wt)
  }

  set.seed(seed)

  dlij <- get_dlij(wt, site_mut, mut_sd_min, dl_sigma)
  f <- get_force(wt, dlij)
  # calculate mutant equilibrium conformation (LRA)
  nzf <- f != 0 # consider only non-zero forces, to make next step faster
  dxyz <-  crossprod(wt$nma$cmat[nzf, ], f[nzf])



  mut <- wt
  mut$graph$lij = wt0$graph$lij + dlij #TODO revise this: mut parameters are w.r.t. w0, not wt...

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

get_dlij <- function(wt, site_mut, mut_sd_min, dl_sigma) {
  graph <- get_graph(wt)

  dlij = rep(0, nrow(get_graph(wt)))

  # pick edges to mutate

  mut_edge <- (graph$i == site_mut | graph$j == site_mut) & (graph$sdij >= mut_sd_min)
  n_mut_edge <- sum(mut_edge)
  dlij[mut_edge] <- rnorm(n_mut_edge, 0, dl_sigma)

  dlij
}


#' Get force resulting from adding dlij to wt
#'
#'
#' @param wt the wild-type protein
#' @param dlij the perturbations to the wt lij parameters
#'
#' @return A force vector of size \code{3 x nsites}
#' @export
#'
#' @examples
#' @family enm mutating functions
get_force <- function(wt, dlij) {

  graph <- get_graph(wt)

  stopifnot(nrow(graph) == length(dlij))

  i <- graph$i
  j <- graph$j
  kij <- graph$kij
  eij <- get_eij(wt)

  fij = -kij * dlij # Force on i in the direction from i to j.


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
