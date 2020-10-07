#' Calculate double-mutational-scan matrix using analytical "fast" method
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


  edges <- seq(nrow(g))
  edges_of_site = list()
  for (p in seq(nsites)) {
    edges_of_site <-  c(edges_of_site, list(edges[g$i == p | g$j == p]))
  }


  for (p in seq(nsites)) {
    for( q in seq(nsites)) {
      drpdrqmin[p, q] <- -sum(abs(dd[edges_of_site[[p]], edges_of_site[[q]]]))
    }
  }

  drpdrqmin

}





# Alternatives ------------------------------------------------------------


#' Calculate site-by-site "fast" response matrix, general
#'
calculate_sisjmin.fast_v1 <- function(wt, amat, mut_dl_sigma, mut_sd_min) {

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



