#' Calculate double-mutational-scanning matrix
#'

sdmrs <- function(wt, nmut, mut_dl_sigma, mut_sd_min, seed, option = "mean_max") {
  # I tested this on a protein with 100 sites and 100 mutations, it's not better than the other option, and it will take more memory...

  cmat <- get_cmat(wt)


  tic()
  fmati <- generate_fmat_dmrs(wt, nmut, mut_dl_sigma, mut_sd_min, 1 * seed)
  fmatj <- generate_fmat_dmrs(wt, nmut, mut_dl_sigma, mut_sd_min, 2 * seed)
  t <-  toc(quiet = T)
  t_fmat <- t$toc - t$tic


  nsites <- get_nsites(wt)

  # structural change due to mutations fmati
  dim(fmati) <- c(3 * nsites, nsites * nmut)
  dri <-  as.matrix(Matrix(cmat) %*% Matrix(fmati, sparse = T))
  dim(dri) <- c(3 * nsites, nsites, nmut)


  # structural change due to mutations fmatj
  dim(fmatj) <- c(3 * nsites, nsites * nmut)
  drj <-  as.matrix(Matrix(cmat) %*% Matrix(fmatj, sparse = T))
  dim(drj) <- c(3 * nsites, nsites, nmut)

  drj2 <- drj^2 %>%
    apply(MARGIN = c(2, 3), FUN = sum)
  dmrs_matrix <- matrix(NA_real_, nrow = nsites, ncol = nsites)

  # compensation matrices

  tic()

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
  t <-  toc(quiet = T)
  t_dridrj <- t$toc - t$tic

  lst(dmrs_matrix, t_fmat, t_dridrj)
}





#' Generate force-matrix for dmrs_matrix calculation
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
      dlmat[, j, m] <- get_delta_lij_dmrs(wt, site_mut = j, mut_sd_min, mut_dl_sigma)
      fmat[, j, m] <- get_force_dmrs(wt, dlmat[, j, m])
    }
  }
  fmat
}




#' get force from delta_lik for dmrs calculations
#'
get_force_dmrs <- function(wt, delta_lij) {

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


get_delta_lij_dmrs <- function(wt, site_mut, mut_sd_min, mut_dl_sigma) {
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
