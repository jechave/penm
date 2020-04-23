#' Get mutant
#'
#' Gets a mutant by introducing a mutation at a random site
#'
#' @param wt The protein \code{prot} to mutate
#' @param wt0 A protein with respect to which mutations are to be calculated
#' @param ideal A protein with ideal protein structure
#' @param param A list containig \code{enm, mut, fit} parameters
#'
#' @return A list of two components, the mutant and the mutated site: \code{lst(mut, site_mut)}
#' @export
#'
#' @examples
#'
#' @family enm mutating functions
get_mutant <- function(wt, wt0, ideal = wt0, param) {
  # pick site to mutate
  site_mut <- sample(wt$site, 1)
  # get mutant
  mut <- get_mutant_site(wt, site_mut, wt0, ideal, param)
  # return mutant and mutated site
  lst(mut = mut, site_mut = site_mut)
}

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
#' @param wt0 A protein with respect to which mutations are to be calculated
#' @param ideal A protein with ideal protein structure
#' @param param A list containig \code{enm, mut, fit} parameters
#'
#' @return
#' @export
#'
#' @examples
#'
#' @family enm mutating functions
#'
get_mutant_site <- function(wt, site_mut , wt0, ideal = wt0, param) {
  # copy wt
  mut <- wt

  # pick edges to mutate
  mut_edge <-
    (mut$enm$graph$i == site_mut |
       mut$enm$graph$j == site_mut) &
    mut$enm$graph$sdij >= param$mut$sd_min

  # mutate edges
  n_mut_edge <- sum(mut_edge)
  dl_sigma <- param$mut$dl_sigma
  dlij <- rnorm(n_mut_edge, 0, dl_sigma)
  #dlij <- sample(c(-dl_max,0,dl_max),size = n_mut_edge, replace = TRUE)
  mut$enm$graph$lij[mut_edge] <-  wt0$enm$graph$lij[mut_edge] + dlij

  # calculate mutant equilibrium conformation (LRA)
  f <- get_force(wt, mut)
  nzf <- f != 0 # consider only non-zero forces
  dxyz <-  crossprod(wt$enm$cmat[nzf, ], f[nzf])
  mut$xyz <- wt$xyz + dxyz

  # recalculate equilibrium distances
  mut$enm$graph$dij <- dij_edge(mut$xyz, mut$enm$graph$i, mut$enm$graph$j)

  # recalculate energies
  mut$energy <- energy(mut, ideal = ideal)

  # revise
  # # recalculate fitness
  # mut$fitness <- fitness(mut,param$fit)
  # mut$log_fit <- log_fitness_prot(mut, param$fit)
  mut
}
