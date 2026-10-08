#' Get a single-point mutant
#'
#' Returns a mutant given a wt and a site to mutate (site_mut)
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation An integer, if 0, return \code{wt} without mutating
#'
#' @return A mutated protein object
#'
#' @details
#' How the mutant is made — the mutational model and its parameters, and the
#' \code{ensemble} — is part of \code{wt}: it was set by [set_enm()], and the
#' mutant inherits it, so that it can be mutated in turn the same way.
#'
#' The mutation is a set of random perturbations of the contacts of
#' \code{site_mut}; there are no amino acids in this model, and no finite set
#' of mutations to draw from. Which perturbations a given mutant gets is fixed
#' by \code{(ensemble, site_mut, mutation)} — see \code{?penm_ensemble} for
#' what that means and when \code{ensemble} may be changed.
#'
#' @export
#'
#' @seealso [set_enm()] to build the `wt` argument and choose how it mutates;
#'   [penm_ensemble] for what `ensemble` means and when to change it;
#'   [delta_structure_by_site], [delta_motion_by_site] and [delta_energy] to
#'   measure the resulting wild-type-vs-mutant differences.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, ensemble = 7)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1)
#'
#' # mutation = 0 returns wt unchanged
#' identical(get_mutant_site(wt, site_mut = 11, mutation = 0), wt)
#'
#' @family enm mutating functions
#'
get_mutant_site <- function(wt, site_mut, mutation = 0) {
  param <- get_enm_param(wt)
  if (is.null(param$mut_model)) {
    stop("wt has no mut_model: it was built by an earlier version of penm. Rebuild it with set_enm().")
  }

  if (param$mut_model == "lfenm") {
    mut <- get_mutant_site_lfenm(wt, site_mut, mutation, param$mut_dl_sigma, param$mut_sd_min, param$ensemble)
    return(mut)
  }

  stop("get_mutant_site does not support mut_model = \"", param$mut_model, "\" yet")

}

#' Get a single-point mutant using lfenm model
#'
#' Returns a mutant given a wt and a site to mutate (site_mut)
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation An integer, if 0, return \code{wt} without mutating
#' @param mut_dl_sigma The standard deviation of a normal distribution from which edge-length perturbation is picked.
#' @param mut_sd_min An integer, only edges with \code{sdij >= mut_sd_min} are mutated
#' @param ensemble An integer naming which realization of the mutational process
#'   the mutant belongs to. With \code{ensemble} fixed, \code{(site_mut, mutation)}
#'   names one specific, reproducible set of contact perturbations. Hold it
#'   constant across a scan or a trajectory; see \code{?penm_ensemble}.
#'
#' @return A mutated protein

#' @noRd
#'
#'
#' @family enm mutating functions
#'
get_mutant_site_lfenm <- function(wt, site_mut, mutation, mut_dl_sigma, mut_sd_min,  ensemble) {

  if (mutation == 0) {
    # if mutation is 0, return wt
    return(wt)
  }

  delta_lij <- with_mut_seed(
    mut_seed(ensemble, site_mut, mutation),
    generate_delta_lij(wt, site_mut, mut_sd_min, mut_dl_sigma)
  )
  f <- calculate_force(wt, delta_lij)
  dxyz <- calculate_dxyz(wt, f)
  mut <- wt
  mut$graph$lij <-  wt$graph$lij + delta_lij #TODO revise this: mut parameters are w.r.t. w0, not wt...
  mut$nodes$xyz <- wt$nodes$xyz + dxyz
  mut$graph$dij <- dij_edge(mut$nodes$xyz, mut$graph$i, mut$graph$j)
  return(mut)
}



#' Perturbations (delta_lij) of contacts of mutated site
#'
#' @noRd
#' @family enm mutating functions
generate_delta_lij <- function(wt, site_mut, mut_sd_min, mut_dl_sigma) {
  graph <- get_graph(wt)

  delta_lij <-  rep(0, nrow(get_graph(wt)))

  # pick edges to mutate

  mut_edge <- (graph$i == site_mut | graph$j == site_mut) & (graph$sdij >= mut_sd_min)
  n_mut_edge <- sum(mut_edge)
  delta_lij[mut_edge] <- rnorm(n_mut_edge, 0, mut_dl_sigma)

  delta_lij
}

#' Calculate structural response of network to applied force
#'
#' \eqn{\delta\mathbf{r} = \mathbf{C}\mathbf[f]}
#'
#' @noRd
#'
#' @family enm mutating functions
calculate_dxyz <- function(wt, f) {
  cmat <- get_cmat(wt)
  nzf <- f != 0 # consider only non-zero forces, to make next step faster
  dxyz <-  crossprod(cmat[nzf, ], f[nzf]) # calculate mutant equilibrium conformation (LRA)
  dxyz
}


#' Get force resulting from adding delta_lij to wt
#'
#'
#' @param wt the wild-type protein
#' @param delta_lij the perturbations to the wt lij parameters
#'
#' @return A force vector of size \code{3 x nsites}
#'
#' @noRd
#'
#'
#' @family enm mutating functions
calculate_force <- function(wt, delta_lij) {

  graph <- get_graph(wt)

  stopifnot(nrow(graph) == length(delta_lij))

  graph$dlij  <- delta_lij

  graph <- graph %>%
    filter(dlij != 0)

  stopifnot(nrow(graph) > 0)



  i <- graph$i
  j <- graph$j
  kij <- graph$kij
  dlij <- graph$dlij

  eij <- get_eij(wt)[delta_lij != 0, ]

  if(nrow(graph) == 1) dim(eij) = c(1, 3)


  fij <-  -kij * dlij # Force on i in the direction from i to j.

  f <- matrix(0, nrow = 3, ncol = get_nsites(wt))

  for (k in seq(nrow(graph))) {
    ik <- i[k]
    jk <- j[k]
    f[, ik] <- f[, ik] + fij[k] * eij[k, ]
    f[, jk] <- f[, jk] - fij[k] * eij[k, ]
  }
  as.vector(f)
}


