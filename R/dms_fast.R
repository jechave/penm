#' Calculate double-mutational-scan matrix using analytical "fast" method
#'
dmrs_analytical <- function(wt, mut_dl_sigma, mut_sd_min, option = "mean_min") {
  stopifnot(option == "mean_min" | option == "min_min")
  cmat <- get_cmat(wt)
  dmrs_analytical_amat(wt, cmat, mut_dl_sigma, mut_sd_min, option)
}




#' Calculate site-by-site "fast" response matrix, general, o
#' minimize over mutations at i and 3yyj
#'
dmrs_analytical_amat <- function(wt, amat, mut_dl_sigma, mut_sd_min, option) {
  stopifnot(option == "mean_min" | option == "min_min")

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

  dd <- t(wd) %*% wd # matrix of size nedges x nedges

  edges <- seq(nrow(g))
  edges_of_site = list()
  cn <- c()
  alpha <- c()
  for (p in seq(nsites)) {
    edges_of_p <- list(edges[g$i == p | g$j == p])
    cn <- c(cn, length(edges_of_p[[1]]))
    alpha <- c(alpha, mut_dl_sigma * sqrt(cn))  #  length of force vector f_p (in edge space)
    edges_of_site <-  c(edges_of_site, edges_of_p)
  }

  dmrs_matrix = matrix(NA, nsites, nsites)
  for (i in seq(nsites)) {
    for( j in seq(nsites)) {
      dd_ij <- dd[edges_of_site[[i]], edges_of_site[[j]]] # edge-edge i-j submatrix
      mij = dd_ij %*% t(dd_ij) # extrama have to do with eigenvectors and values of matrix m = a * tr(a)
      if (option == "mean_min")
        dmrs_matrix[i, j] <-  - alpha[i] * alpha[j] * sqrt(sum(diag(mij)) / cn[i]) # sqrt(mean(min(dri.drj)^2))
      if (option == "min_min")
        dmrs_matrix[i, j] <- - alpha[i] * alpha[j] * sqrt(eigen(mij)$values[1]) # minimum dri.drj
    }
  }

  dmrs_matrix

}

