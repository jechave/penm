#' Get a single-point mutant
#'
#' Returns a mutant given a wt and a site to mutate (site_mut)
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation The allele \code{site_mut} is given: an integer from 0 to
#'   \code{n_alleles - 1} (see [set_enm()]). Allele 0 is the pdb's residue. If
#'   the site already has this allele, \code{wt} is returned unchanged.
#'
#' @return A mutated protein object
#'
#' @details
#' How the mutant is made — the mutational model and its parameters, and the
#' \code{ensemble} — is part of \code{wt}: it was set by [set_enm()], and the
#' mutant inherits it, so that it can be mutated in turn the same way.
#'
#' Each site carries an allele. There are no amino acids in this model: what an
#' allele does is a set of random perturbations of the equilibrium lengths of
#' the site's contacts, fixed by \code{(ensemble, site_mut, mutation)} — see
#' \code{?penm_ensemble} for what that means and when \code{ensemble} may be
#' changed. The equilibrium lengths depend only on the alleles, so mutating a
#' site back to an earlier allele restores its earlier lengths.
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
#' # the allele a site already has returns the protein unchanged
#' identical(get_mutant_site(wt, site_mut = 11, mutation = 0), wt)
#'
#' # mutating back restores the wild type's equilibrium lengths
#' back <- get_mutant_site(mut, site_mut = 11, mutation = 0)
#' identical(back$graph$lij, wt$graph$lij)
#'
#' @family enm mutating functions
#'
get_mutant_site <- function(wt, site_mut, mutation = 0) {
  param <- get_enm_param(wt)
  if (is.null(param$mut_model) || is.null(wt$nodes$sequence)) {
    stop("wt has no mut_model or no sequence: it was built by an earlier version of penm. Rebuild it with set_enm().")
  }
  nsites <- get_nsites(wt)
  if (!(site_mut %in% seq_len(nsites))) {
    stop("site_mut must be one of 1..", nsites)
  }
  if (!(mutation %in% 0:(param$n_alleles - 1))) {
    stop("mutation must be an allele, one of 0..", param$n_alleles - 1)
  }
  if (mutation == wt$nodes$sequence[site_mut]) {
    # the protein with this allele at site_mut is wt itself
    return(wt)
  }

  if (param$mut_model == "lfenm") {
    mut <- get_mutant_site_lfenm(wt, site_mut, mutation)
    return(mut)
  }

  stop("get_mutant_site does not support mut_model = \"", param$mut_model, "\" yet")

}

#' Get a single-point mutant using lfenm model
#'
#' Gives `site_mut` the allele `mutation`, recomputes the equilibrium lengths,
#' and moves the structure by the linear response to the force their change
#' exerts. kmat, the normal modes and the edge directions are those of the
#' protein `set_enm()` built, so the mutant's structure is that protein's plus
#' the response to `lij - l0ij`, whatever path of mutations led to it.
#'
#' @param wt The protein \code{prot} to mutate
#' @param site_mut The site to mutate (not the pdb_site, but sequential)
#' @param mutation The new allele, different from the site's current one
#'
#' @return A mutated protein

#' @noRd
#'
#'
#' @family enm mutating functions
#'
get_mutant_site_lfenm <- function(wt, site_mut, mutation) {
  mut <- wt
  mut$nodes$sequence[site_mut] <- as.integer(mutation)
  mut$graph$lij <- calculate_lij(mut)
  delta_lij <- mut$graph$lij - wt$graph$lij
  f <- calculate_force(wt, delta_lij)
  dxyz <- calculate_dxyz(wt, f)
  mut$nodes$xyz <- wt$nodes$xyz + dxyz
  mut$graph$dij <- dij_edge(mut$nodes$xyz, mut$graph$i, mut$graph$j)
  return(mut)
}


#' Equilibrium lengths implied by the sequence
#'
#' \eqn{l_{ij} = l^0_{ij} + \delta(i, s_i)_{ij} + \delta(j, s_j)_{ij}} for every
#' edge, with \eqn{\delta} given by [allele_delta_lij()] for edges with
#' `sdij >= mut_sd_min`, and 0 for the others. Computed from `l0ij` and the
#' sequence only, so equal sequences give identical lengths. Used by both
#' mutational models.
#'
#' @param prot a `prot`
#'
#' @returns `lij`, one per edge
#'
#' @noRd
#' @family enm mutating functions
#'
calculate_lij <- function(prot) {
  graph <- prot$graph
  sequence <- prot$nodes$sequence
  perturbed <- graph$sdij >= prot$param$mut_sd_min
  delta_from_i <- numeric(nrow(graph))
  delta_from_j <- numeric(nrow(graph))

  for (site in which(sequence != 0)) {
    delta_lij <- allele_delta_lij(prot, site, sequence[site])
    site_is_i <- perturbed & graph$i == site
    site_is_j <- perturbed & graph$j == site
    delta_from_i[site_is_i] <- delta_lij[graph$j[site_is_i]]
    delta_from_j[site_is_j] <- delta_lij[graph$i[site_is_j]]
  }

  # always added in this order: floating-point addition is not associative
  (graph$l0ij + delta_from_i) + delta_from_j
}


#' The change an allele makes to the edges of its site
#'
#' A vector of `nsites` normal draws with sd `mut_dl_sigma`, seeded by
#' `(ensemble, site, allele)` without disturbing the caller's RNG. Element `k`
#' is the change to the edge between `site` and site `k` (element `site` itself
#' is not used). Allele 0 changes nothing.
#'
#' @param prot a `prot`
#' @param site sequential site index
#' @param allele an allele in `0 .. n_alleles - 1`
#'
#' @returns a vector of length `nsites`
#'
#' @noRd
#' @family enm mutating functions
#'
allele_delta_lij <- function(prot, site, allele) {
  nsites <- get_nsites(prot)
  if (allele == 0) return(numeric(nsites))
  seed <- mut_seed(prot$param$ensemble, site, allele)
  with_mut_seed(seed, stats::rnorm(nsites, mean = 0, sd = prot$param$mut_dl_sigma))
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
  as.vector(dxyz) # crossprod gives a 3N x 1 matrix; xyz is a vector
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


