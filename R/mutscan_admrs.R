#' Calculate a double-mutational-scan matrix analytically
#'
#' Returns a compensation matrix: element (i,j) measures the degree of compensation of structural deformations produced by pairs of mutations at sites i and j.
#' It uses analytical methods (closed formulas).
#' Two measures are implemented:
#' "mean_max" (default), the structural compensation maximized over mutations at j and averaged over mutations at i;
#' "max_max" is the structural compensation maximized over mutations at i and j.
#'
#' For details see \doi{10.7717/peerj.11330}
#'
#' @param wt is the (wild-type) protein to mutate (an object obtained using \code{set_enm})
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
#' dmat <- admrs(wt, mut_dl_sigma = 0.3, mut_sd_min = 1, option = "max_max", response = "dr2")
#' }
#'
#' @family mutscan functions
#'
admrs <- function(wt, mut_dl_sigma, mut_sd_min, option = "mean_max", response = "dr2") {
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
  admrs_amat(wt, amat, mut_dl_sigma, mut_sd_min, option)
}




#' Calculate site-by-site analytic response matrix
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
admrs_amat <- function(wt, amat, mut_dl_sigma, mut_sd_min, option) {
  stopifnot(option == "mean_max" | option == "max_max")

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
      mij = dd_ij %*% t(dd_ij) # extrema have to do with eigenvectors and values of matrix m = a * tr(a)
      if (option == "mean_max")
        dmrs_matrix[i, j] <-   alpha[i] * alpha[j] * sqrt(sum(diag(mij)) / cn[i]) # sqrt(mean(max(dri.drj)^2))
      if (option == "max_max")
        dmrs_matrix[i, j] <-  alpha[i] * alpha[j] * sqrt(eigen(mij)$values[1]) #  maximum of dri.drj
    }
  }

  dmrs_matrix

}

