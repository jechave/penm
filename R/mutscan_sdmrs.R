#' Calculate a double-mutational-scan matrix numerically (simulation-based)
#'
#' Returns a compensation matrix: element (i,j) measures the degree of compensation of structural deformations produced by pairs of mutations at sites i and j.
#' It uses a simulation method (calculates responses for various instances of forces, then calculates means or maxima)
#' Two measures are implemented:
#' "mean_max" (default), the structural compensation maximized over mutations at j and averaged over mutations at i;
#' "max_max" is the structural compensation maximized over mutations at i and j.
#'
#' For details see \doi{10.7717/peerj.11330}
#'
#' @param wt is the (wild-type) protein to mutate (an object obtained using \code{set_enm})
#' @param nmut is the number of mutations per site to simulate
#' @param mut_dl_sigma is the standard deviation of a normal distribution from which edge-length perturbations are picked (LFENM model).
#' @param mut_sd_min is integer sequence-distance cutoff, only edges with \code{sdij >= mut_sd_min} are mutated
#' @param option is either "mean_max" (default) or "max_max", depending on which compensation measure is desired.
#' @param response is the response desired, which maybe either "dr2", "de2", or "df2"
#'
#' @return A compensation matrix, rows are initially mutated site, j is compensation site
#' @export
#'
#' @examples
#' \dontrun{
#' pdb <- bio3d::read.pdb("2acy")
#' wt <- set_enm(pdb, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#' dmat <- sdmrs(wt, nmut = 10, mut_dl_sigma = 0.3, mut_sd_min = 1, option = "max_max", response = "dr2")
#' }
#'
#' @family mutscan functions
#'
sdmrs <- function(wt, nmut, mut_dl_sigma, mut_sd_min,  option = "mean_max", response = "dr2", seed = 1024) {
  stopifnot(option == "mean_max" | option == "max_max")

  if (response == "dr2") {
    amat <- get_cmat(wt)
  } else if (response == "de2") {
    amat <- get_cmat_sqrt(wt)
  } else if (response == "df2") {
    amat <- diag(3 * get_nsites(wt))
  } else {
    stop("Unknown value of response, stop admrs")
  }
  result <- sdmrs_amat(wt, amat, nmut, mut_dl_sigma, mut_sd_min, option, seed)
  result
}



#' Calculate site-by-site response matrix
#'
#' This is a general function for a response vector of the form \code{amat * f}.
#' If \code{amat = cmat}, then it's the structural deformation,
#' for \code{amat = identity}, then it's the force vector,
#' for \code{amat = cmat^{1/2}}, it's a deformation-energy vector.
#'
#' @noRd
#'
#' @family mutscan functions
#'
sdmrs_amat <- function(wt, amat, nmut, mut_dl_sigma, mut_sd_min,  option, seed) {

  stopifnot(option == "mean_max" | option == "max_max")

  fmati <- generate_fmat_dmrs(wt, nmut, mut_dl_sigma, mut_sd_min, 1 * seed)
  fmatj <- generate_fmat_dmrs(wt, nmut, mut_dl_sigma, mut_sd_min, 2 * seed)


  nsites <- get_nsites(wt)

  # response vectors due to mutations fmati (if amat = cmat, then response is structural change)
  dim(fmati) <- c(3 * nsites, nsites * nmut)
  dri <-  as.matrix(Matrix::Matrix(amat) %*% Matrix::Matrix(fmati, sparse = T))
  dim(dri) <- c(3 * nsites, nsites, nmut)


  # response vectors due to mutations fmati (if amat = cmat, then response is structural change)
  dim(fmatj) <- c(3 * nsites, nsites * nmut)
  drj <-  as.matrix(Matrix::Matrix(amat) %*% Matrix::Matrix(fmatj, sparse = T))
  dim(drj) <- c(3 * nsites, nsites, nmut)

  drj2 <- drj^2 %>%
    apply(MARGIN = c(2, 3), FUN = sum)
  dmrs_matrix <- matrix(NA_real_, nrow = nsites, ncol = nsites)

  # compensation matrices

  for (i in seq(nsites)) {
    for (j in seq(nsites)) {
      mi <- dri[ , i, ]
      mj <- drj[ , j, ]
      mimj <- t(mi) %*% mj
      if (option == "mean_max")
        dmrs_matrix[i,j] <- sqrt(mean(matrixStats::rowMaxs(mimj^2)))
      if (option == "max_max")
        dmrs_matrix[i,j] <- max(mimj^2)
    }
  }

  dmrs_matrix
}





#' Generate force-matrix for dmrs_matrix calculation
#'
#' @noRd
#'
#' @family mutscan functions
#'
generate_fmat_dmrs <- function(wt, nmut, mut_dl_sigma, mut_sd_min,  seed) {
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
      dlmat[, j, m] <- generate_delta_lij_dmrs(wt, site_mut = j, mut_sd_min, mut_dl_sigma)
      fmat[, j, m] <- calculate_force_dmrs(wt, dlmat[, j, m])
    }
  }
  fmat
}




#' get force from delta_lik for dmrs calculations
#'
#' @noRd
#'
#' @family mutscan functions
#'
calculate_force_dmrs <- function(wt, delta_lij) {

  graph <- get_graph(wt)

  stopifnot(nrow(graph) == length(delta_lij))

  graph$dlij  <- delta_lij

  mutated_edge <- !near(delta_lij, 0)
  graph <- graph[mutated_edge, ]


  i <- graph$i
  j <- graph$j
  kij <- graph$kij
  dlij <- graph$dlij

  fij <-   dlij  # here dlij is a constant force, in LFENM it should be kij * dlij. Identical for ANM.

  eij <- get_eij(wt)[mutated_edge, ]


  f <- matrix(0., nrow = 3, ncol = get_nsites(wt))

  for (k in seq(nrow(graph))) {
    ik <- i[k]
    jk <- j[k]
    f[, ik] <- f[, ik] + fij[k] * eij[k, ]
    f[, jk] <- f[, jk] - fij[k] * eij[k, ]
  }
  as.vector(f)
}


#' get spring-length  perturbations dmrs calculations
#'
#' @noRd
#'
#' @family mutscan functions
#'
generate_delta_lij_dmrs <- function(wt, site_mut, mut_sd_min, mut_dl_sigma) {
  graph <- get_graph(wt)

  delta_lij <-  rep(0, nrow(get_graph(wt)))

  # pick edges to mutate

  mut_edge <- (graph$i == site_mut | graph$j == site_mut) & (graph$sdij >= mut_sd_min)
  n_mut_edge <- sum(mut_edge)

  # dl <- runif(n_mut_edge, -sqrt(3) * mut_dl_sigma, sqrt(3) * mut_dl_sigma)
  dl <- rnorm(n_mut_edge, 0, mut_dl_sigma)
  dl <- dl * mut_dl_sigma / sqrt(mean(dl^2)) # constraint: |dl|^2 = mut_dl_sigma^2 * N_edges

  delta_lij[mut_edge] <- dl

  delta_lij
}



