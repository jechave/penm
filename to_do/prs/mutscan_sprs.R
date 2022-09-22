#' Calculate Bahar's PRS matrix, simulation-based
#'
#' @param cmat is a covariance matrix
#' @param n_randforces is the number of random forces to apply to each site
#' @param normalize is a logical valuable of whether to normalize the PRS matrix (default is F)
#'
#' @returns
#'
#' @export
#'
sprs <- function(covar_mat, n_randforces, normalize = F){
  a <- -0.5
  b <- 0.5
  n_res <- (nrow(covar_mat))/3
  PRS_matrix <- matrix(0, n_res, n_res)
  for (j in 1:n_res){
    rand_forces <- a + (b-a) * matrix(runif(3 * n_randforces), 3, n_randforces)
    colnorm <- sqrt(colSums(rand_forces^2)) # Se generan las fuerzas al azar
    rand_forces <- rand_forces/(rbind(rbind(colnorm, colnorm), colnorm)) # Se normalizan las fuerzas
    force_mat <- matrix (0, 3 * n_res, n_randforces)
    force_mat[3 * j - 2,] <- rand_forces[1,]
    force_mat[3 * j - 1,] <- rand_forces[2,]
    force_mat[3 * j,] <- rand_forces[3,]
    deltaR <- (covar_mat %*% force_mat)^2
    deltaRx <- deltaR[seq(1, nrow(deltaR), by = 3),]
    deltaRy <- deltaR[seq(2, nrow(deltaR), by = 3),]
    deltaRz <- deltaR[seq(3, nrow(deltaR), by = 3),]
    PRS_matrix[, j] <- rowSums(deltaRx + deltaRy + deltaRz)/n_randforces
  }


  if (normalize) {
    diagonal <- diag(1 /diag(PRS_matrix))
    normal <- diagonal %*% PRS_matrix # Se normaliza dividiendo la matriz por el elemneto de la diagonal en la misma columna
    normal_PRS <- normal - diag(n_res) # Los elementos de la diagonal se hacen 0
    result <- normal
  } else {
    result <- PRS_matrix
  }
  result
}




