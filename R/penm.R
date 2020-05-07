#' Get a single-point mutant
#'
#' Returns a mutant given a wt and a site to mutate (site_mut)
#'     wt0 is the protein to use to calculate perturbations of current site
#'     (if wt0 = wt, it is the current wild-type, but it can be a fixed initial state)
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation An integer, if 0, return \code{wt} without mutating
#' @param seed An integer, the seed for set.seed before picking perturbations
#' @param wt0 A protein with respect to which mutations are to be calculated
#' @param mut_sd_min An integer, only edges with \code{sdij > mut_sd_min} are mutated
#' @param dl_sigma The standard deviation of a normal distribution from which edge-length perturbation is picked.
#' @param model The enm model type
#' @param d_max The enm distance cut-off to define contacts
#' @param frustrated Whether to include frustrations or not
#' @param update_enm Whether to update enm or not
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
                            mut_sd_min = 2, dl_sigma = .3, update_enm = FALSE,
                            model = "ming_wall", d_max = 10.5, frustrated = FALSE) {

  if (mutation == 0) {
    # if mutation is 0, return wt
    return(wt)
  }

  set.seed(seed)

  mut <- wt

  # pick edges to mutate
  mut_edge <-
    (mut$enm$graph$i == site_mut |
       mut$enm$graph$j == site_mut) &
    mut$enm$graph$sdij >= mut_sd_min

  # mutate edges
  n_mut_edge <- sum(mut_edge)
  dlij <- rnorm(n_mut_edge, 0, dl_sigma)
  #dlij <- sample(c(-dl_max,0,dl_max),size = n_mut_edge, replace = TRUE)
  mut$enm$graph$lij[mut_edge] <-  wt0$enm$graph$lij[mut_edge] + dlij # add perturbation to spring parameters

  # calculate mutant equilibrium conformation (LRA)
  f <- get_force(wt, mut)
  nzf <- f != 0 # consider only non-zero forces, to make next step faster
  dxyz <-  crossprod(wt$enm$cmat[nzf, ], f[nzf])
  mut$xyz <- wt$xyz + dxyz

  # recalculate enm

  if (update_enm) { # recalculate enm
    # recalculate whole mut$enm (graph, kmat, etc. etc.)
    mut <- mutate_enm(mut, model, d_max, frustrated)
  } else {
  # recalculate only equilibrium distances
    mut$enm$graph$dij <- dij_edge(mut$xyz, mut$enm$graph$i, mut$enm$graph$j)
  }

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





#' mutate enm following change in protein structure and lij parameters
#'
mutate_enm <- function(prot, model, d_max, frustrated) {
  # it re-calculates what needs to be recalculated due to change in graph


  prot <- mutate_graph(prot, model, d_max) # recalculate mutant's graph
  prot <- mutate_eij(prot) # recalculate eij versors, following change in graph
  prot$enm$kmat <- set_enm_kmat(prot$enm$graph, prot$enm$eij, prot$nsites, frustrated) # recalculate kmat

  # recalculate normal modes etc.
  nma <- enm_nma(prot$enm$kmat) #returns mode, evalue, cmat, umat

  prot$enm$mode <- nma$mode
  prot$enm$evalue <- nma$evalue
  prot$enm$cmat <- nma$cmat
  prot$enm$umat <- nma$umat

  prot
}

#' mutate graph following change in structure
#'
mutate_graph <- function(mut, model, d_max) {
  g1 <- mut$enm$graph # the mut graph with wt contacts but mut lij parameters...

  # the "self-consistent" graph: add/rm edges according to new xyz
  g2 <- set_enm_graph(mut$xyz, mut$pdb_site, model = model,  d_max = d_max)

  # the "frustrated" graph: keep lij for edges that haven't changed
  g2[g2$edge %in% g1$edge, "lij"] <- g1[g1$edge %in% g2$edge, "lij"]

  mut$enm$graph <- g2
  mut
}

#' mutate eij vectors, following change of structure and graph
mutate_eij <- function(mut) {
  mut$enm$eij <- set_enm_eij(mut$xyz, mut$enm$graph$i, mut$enm$graph$j)
  mut
}

#' Depreciated, use mutate_enm instead
#'
enm_update <- function(...) {
  stop("enm_update has been renamed to mutate_enm")
}



