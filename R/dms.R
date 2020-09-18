#' Calculate double-mutational-scanning matrix
#'
dmsmat_dms.new <- function(wt, nmut, mut_dl_sigma, mut_sd_min, seed) {
  cmat <- get_cmat(wt)
  fmati <- generate_fmat_dms(wt, nmut, mut_dl_sigma, mut_sd_min, 1 * seed)
  fmatj <- generate_fmat_dms(wt, nmut, mut_dl_sigma, mut_sd_min, 2 * seed)

  nsites <- get_nsites(wt)

  # structural change due to mutations fmati
  dim(fmati) <- c(3 * nsites, nsites * nmut)
  dri <-  as.matrix(Matrix(cmat) %*% Matrix(fmati, sparse = T))
  dim(dri) <- c(3 * nsites, nsites, nmut)


  # structural change due to mutations fmatj
  dim(fmatj) <- c(3 * nsites, nsites * nmut)
  drj <-  as.matrix(Matrix(cmat) %*% Matrix(fmatj, sparse = T))
  dim(drj) <- c(3 * nsites, nsites, nmut)

  dmsmat <- matrix(NA, nrow = nsites, ncol = nsites)
  for (i in seq(nsites)) {
    for (j in seq(nsites)) {
      dridrj <- t(dri[, i, ]) %*% drj[, j, ]
      dmsmat[i,j] <- -max(abs(dridrj))
    }
  }
  dmsmat

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
      dlmat[, j, m] = sign(dlmat[, j, m])
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

  graph <- graph %>%
    filter(dlij != 0)

  i <- graph$i
  j <- graph$j
  kij <- graph$kij
  dlij <- graph$dlij

  eij <- get_eij(wt)[delta_lij != 0, ]

  fij <-   sign(dlij) # Force on i in the direction from i to j.

  f <- matrix(0, nrow = 3, ncol = get_nsites(wt))

  for (k in seq(nrow(graph))) {
    ik <- i[k]
    jk <- j[k]
    f[, ik] <- f[, ik] + fij[k] * eij[k, ]
    f[, jk] <- f[, jk] - fij[k] * eij[k, ]
  }
  as.vector(f)
}






#' Calculate double-mutational-scan matrix using analytical "fast" method
#'
#' Calculate site-by-site structure response matrix dr2(n, j) using method "fast"
#'

dmsmat_dms.fast <- function(wt, mut_dl_sigma, mut_sd_min) {
  cmat <- get_cmat(wt)
  calculate_sisjmin.fast(wt, cmat, mut_dl_sigma, mut_sd_min)
}

' Calculate site-by-site "fast" response matrix, general
#'

calculate_sisjmin.fast <- function(wt, amat, mut_dl_sigma, mut_sd_min) {

  g <- get_graph(wt)
  eij  <- get_eij(wt)

  nsites <- get_nsites(wt)
  dim(amat) <- c(nrow(amat), 3, nsites)

  # response to edge perturbations
  nedges <- nrow(g)
  dmat = matrix(NA, 3 * nsites, nedges)
  for (edge in seq(nedges)) {
    i = g$i[edge]
    j = g$j[edge]
    dmat[, edge] <- (amat[,,j] - amat[,,i]) %*% eij[edge, ]
  }

  weight <-  ifelse(g$sdij >= mut_sd_min, 1, 0)
  wd <- t(t(dmat) * weight) # note that * mutliplies by row

  dd <- t(wd) %*% wd

  drpdrqmin = matrix(NA, nsites, nsites)


  for (p in seq(nsites)) {
    edges_of_p <- g$i == p | g$j == p
    for( q in seq(nsites)) {
      edges_of_q <- g$i == q | g$j == q
      drpdrqmin[p, q] <- -sum(abs(dd[edges_of_p, edges_of_q]))
    }
  }

  drpdrqmin

}





