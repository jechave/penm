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
#' @param seed An integer, the seed for set.seed before picking perturbations
#' @param wt0 A protein with respect to which mutations are to be calculated
#' @param ideal A protein with ideal protein structure
#' @param sd_min An integer, only edges with \code{sdij > sd_min} are mutated
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
                            wt0 = wt, ideal = wt0,
                            sd_min = 2, dl_sigma = .3, update_enm = FALSE,
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


    # recalculate eij versors
    mut$enm$eij <- eij_edge(mut$xyz, mut$enm$graph$i, mut$enm$graph$j)

    # recalculate kmat
    mut$enm$kmat <-  kmat_graph(mut$enm$graph, mut$enm$eij,
                                nsites = length(mut$pdb_site),
                                frustrated = frustrated)
    # do nma and update nma info in mut
    nma <- enm_nma(mut$enm$kmat)
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
#' @param frustrated A logical indicating whether to add frustration
#'
#' @return A protein object with updated \code{enm} component
#' @export
#'
#' @examples
#'
#' @family enm mutating functions
#' @export
enm_update <- function(prot, frustrated = FALSE) {
  # it re-calculates what needs to be recalculated due to change in graph
  stopifnot(!is.null(prot$ind_active)) # stop if no active-diste info
  enm <- prot$enm
  enm$kmat <- kmat_graph(prot$enm$graph, prot$enm$eij, prot$nsites, frustrated)
  nma <- enm_nma(enm$kmat) #returns mode, evalue, cmat, umat
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
