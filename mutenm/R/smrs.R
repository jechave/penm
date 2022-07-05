#' Calcualte dr2ij response matrix, smrs method
#'
smrs <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed) {

  # calculate fmat
  tic()
  perturbations <- generate_perturbations(wt, nmut_per_site, mut_dl_sigma, mut_sd_min, seed)
  fmat <- perturbations$fmat
  t <- toc(quiet = T)
  t_fmat <- t$toc - t$tic

  # calculate dr2ij
  tic()

  dr2ij <- calculate_dr2ij_smrs(wt, fmat)

  t <- toc(quiet = T)
  t_dr2ij <- t$toc - t$tic

  lst(dr2ij,  t_fmat, t_dr2ij)
}


# Generate perturbations (forces) -----------------------------------------

#' Generate matrices of perturbations (dlmat) and forces (fmat)
#'
#' dlmat(edge, j, m) is the perturbation dl at edge e due to mutation m at site j
#' fmat(i, j, m) is force at site i for mutation m = c(1,...nmut), at site j
#'
generate_perturbations <- function(wt, nmut, mut_dl_sigma, mut_sd_min,  seed) {

  mutation = seq(nmut)
  nsites <- get_nsites(wt)
  nedges <- nrow(get_graph(wt))

  dlmat <- matrix(NA, nedges, nsites * nmut)
  dim(dlmat) <- c(nedges, nsites, nmut)

  fmat <- matrix(NA, 3 * nsites, nsites * nmut)
  dim(fmat) <- c(3 * nsites, nsites, nmut)

  for(j in seq(nsites)) {
    for (m in mutation) {
      set.seed(seed + j * m)
      dlmat[, j, m] <- get_delta_lij(wt, site_mut = j, mut_sd_min, mut_dl_sigma)
      fmat[, j, m] <- get_force(wt, dlmat[, j, m])
    }
  }


  lst(dlmat, fmat)

}



#' Calculate structure response matrix dr2(i, j) using method smrs
#'

calculate_dr2ij_smrs <- function(wt, fmat) {
  calculate_s2ij(fmat, get_cmat(wt)) %>%
    matrix_to_tibble() %>%
    rename(dr2ij = mij)
}



#' Calculate response matrices given fmat
#'

calculate_s2ij<- function(fmat, amat) {
  nsites <- dim(fmat)[2]
  nmut <- dim(fmat)[3]

  dim(fmat) <- c(3 * nsites, nsites * nmut)

  smat <- as.matrix(Matrix(amat) %*% Matrix(fmat, sparse = T))

  dim(smat) = c(3, nsites, nsites, nmut)

  s2ij <- smat^2 %>%
    apply(MARGIN = c(2, 3, 4), FUN = sum)  %>%
    apply(MARGIN = c(1, 2), FUN = mean)

  s2ij
}



#'  Get spring length perturbations (to calculate forces afterwars)
#'
get_delta_lij <- function(wt, site_mut, mut_sd_min, mut_dl_sigma) {
  graph <- get_graph(wt)

  delta_lij <-  rep(0, nrow(get_graph(wt)))

  # pick edges to mutate

  mut_edge <- (graph$i == site_mut | graph$j == site_mut) & (graph$sdij >= mut_sd_min)
  n_mut_edge <- sum(mut_edge)
  delta_lij[mut_edge] <- rnorm(n_mut_edge, 0, mut_dl_sigma)

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

  graph$dlij  <- delta_lij

  graph <- graph %>%
    filter(dlij != 0)

  i <- graph$i
  j <- graph$j
  kij <- graph$kij
  dlij <- graph$dlij

  eij <- get_eij(wt)[delta_lij != 0, ]

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
