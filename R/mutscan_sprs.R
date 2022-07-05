#' Calculate Bahar's PRS matrix numerically
#'
sprs <- function(covar_mat, n_randforces){
  a <- -0.5
  b <- 0.5
  n_res <- (nrow(covar_mat))/3
  PRS_matrix <- matrix(0, n_res, n_res)
  for (i in 1:n_res){
    rand_forces <- a + (b-a) * matrix(runif(3 * n_randforces), 3, n_randforces)
    colnorm <- sqrt(colSums(rand_forces^2)) # Se generan las fuerzas al azar
    rand_forces <- rand_forces/(rbind(rbind(colnorm, colnorm), colnorm)) # Se normalizan las fuerzas
    force_mat <- matrix (0, 3 * n_res, n_randforces)
    force_mat[3 * i - 2,] <- rand_forces[1,]
    force_mat[3 * i - 1,] <- rand_forces[2,]
    force_mat[3 * i,] <- rand_forces[3,]
    deltaR <- (covar_mat %*% force_mat)^2
    deltaRx <- deltaR[seq(1, nrow(deltaR), by = 3),]
    deltaRy <- deltaR[seq(2, nrow(deltaR), by = 3),]
    deltaRz <- deltaR[seq(3, nrow(deltaR), by = 3),]
    PRS_matrix[, i] <- rowSums(deltaRx + deltaRy + deltaRz)/n_randforces
  }

  diagonal <- diag(1 /diag(PRS_matrix))
  normal <- diagonal %*% PRS_matrix # Se normaliza dividiendo la matriz por el elemneto de la diagonal en la misma columna
  normal_PRS <- normal - diag(n_res) # Los elementos de la diagonal se hacen 0

  #write.table(normal_PRS, "Normal_PRS.dat", sep = "\t", col.names = FALSE, row.names = FALSE)
  #write.table(PRS_matrix, "PRS.dat", sep = "\t", col.names = FALSE, row.names = FALSE)
  return(list(PRS_matrix = PRS_matrix, normal_PRS = normal_PRS))
}

# Perfiles de sensibilidad e influencia
perfiles <- function(m) {
  sensibilidad <- sqrt((colMeans(m)))
  influencia <- sqrt((rowMeans(m)))
  return(list(sensibilidad = sensibilidad, influencia = influencia))
}


