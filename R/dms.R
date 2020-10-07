
#' Calculate double-mutational-scanning matrix
#'
dmsmat_dms.new <- function(wt, nmut, mut_dl_sigma, mut_sd_min, seed) {
  # I tested this on a protein with 100 sites and 100 mutations, it's not better than the other option, and it will take more memory...

  cmat <- get_cmat(wt)


  tic()
  fmati <- generate_fmat_dms(wt, nmut, mut_dl_sigma, mut_sd_min, 1 * seed)
  fmatj <- generate_fmat_dms(wt, nmut, mut_dl_sigma, mut_sd_min, 2 * seed)
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

  tic()
  dmsmat <- matrix(NA_real_, nrow = nsites, ncol = nsites)
  for (i in seq(nsites)) {
    for (j in seq(nsites)) {
      mi <- dri[ , i, ]
      mj <- drj[ , j, ]
      dridrj <- t(mi) %*% mj
      # dridrj <- t(dri[, i, ]) %*% drj[, j, ] # debug
      dmsmat[i,j] <- -max(abs(dridrj)) # this way I consider fj and -fj and fi and -fi...
    }
  }

  t <-  toc(quiet = T)
  t_dridrj <- t$toc - t$tic


  lst(dmsmat, t_fmat, t_dridrj)
}





#' Generate force-matrix for dmsmat calculation
#'
generate_fmat_dms <- function(wt, nmut, mut_dl_sigma, mut_sd_min,  seed) {
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
      # dlmat[, j, m] = sign(dlmat[, j, m])
      fmat[, j, m] <- get_force_dms(wt, dlmat[, j, m])
    }
  }
  fmat
}




#' get force from delta_lik for dms calculations
#'
get_force_dms <- function(wt, delta_lij) {

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






