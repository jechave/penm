library(here)
## ============================================================================
## WITH vs WITHOUT FRUSTRATION
##
## A network is FRUSTRATED when l_ij != d^e_ij: springs are stressed at the
## equilibrium conformation. Two questions, kept separate:
##   (A) does frustration change the HESSIAN (hence dynamics)?  [the flag]
##   (B) do MUTATIONS create frustration?                       [the model]
## ============================================================================
source(here("R", "penm_catalog.R"))
library(ggplot2)

## penm blocks frustrated=TRUE in set_enm(), so build the Hessian directly.
## Convention (verified earlier): g = l/d - 1 with bracket (ee' - I).
## NOTE: get_eij(prot) is NOT updated by get_mutant_site_lfenm() -- LFENM never
## rebuilds K, so penm has no reason to refresh it. Recompute from xyz, or the
## Hessian mixes wild-type directions with mutant distances.
hess <- function(prot, frustrated) {
  g <- get_graph(prot); N <- get_nsites(prot)
  xyz <- matrix(as.vector(get_xyz(prot)), nrow = 3)
  e <- calculate_enm_eij(xyz, g$i, g$j)
  g$dij <- dij_edge(xyz, g$i, g$j)
  K <- matrix(0, 3*N, 3*N); I3 <- diag(3)
  for (q in seq_len(nrow(g))) {
    i <- g$i[q]; j <- g$j[q]; ee <- tcrossprod(e[q,])
    gg <- if (frustrated) g$lij[q]/g$dij[q] - 1 else 0
    Kij <- -g$kij[q]*(ee + gg*(ee - I3))
    ii <- (3*i-2):(3*i); jj <- (3*j-2):(3*j)
    K[ii,jj] <- Kij; K[jj,ii] <- t(Kij)
  }
  for (i in seq_len(N)) { ii <- (3*i-2):(3*i); b <- matrix(0,3,3)
    for (j in seq_len(N)) if (j!=i) b <- b - K[ii,(3*j-2):(3*j)]; K[ii,ii] <- b }
  K
}
spec <- function(K, tol=1e-8) { e <- eigen(K, symmetric=TRUE)
  keep <- e$values > tol*max(abs(e$values))
  list(val=e$values[keep], vec=e$vectors[,keep,drop=FALSE], nzero=sum(!keep),
       minval=min(e$values), raw=e$values) }

## A saddle test that CAN fail. spec() filters by tolerance, so min(spec$val) is
## structurally positive and tests nothing. Instead: sort the RAW spectrum, drop
## the 6 lowest by RANK (the rigid-body modes), and report the 7th. If the
## network is a saddle that value is negative.
seventh_eigen <- function(K) sort(eigen(K, symmetric=TRUE)$values)[7]
ts_of <- function(v, beta=beta_boltzmann()) sum(0.5/beta*(log(2*pi/(beta*v))+1))

## Relax a prot to the TRUE minimum of its own (frustrated) potential.
## The LFENM structure is the linear-response estimate and is NOT stationary:
## a Hessian evaluated there is not singular on rotations, which produces three
## spurious small eigenvalues. Always relax before any spectral quantity.
enm_force <- function(prot) {
  g <- get_graph(prot); N <- get_nsites(prot); R <- matrix(as.vector(get_xyz(prot)), nrow = 3)
  D <- R[, g$j, drop=FALSE] - R[, g$i, drop=FALSE]
  d <- sqrt(colSums(D^2)); E <- sweep(D, 2, d, "/")
  a <- g$kij * (d - g$lij); f <- matrix(0, 3, N)
  for (q in seq_len(nrow(g))) { f[, g$i[q]] <- f[, g$i[q]] + a[q]*E[,q]
                                f[, g$j[q]] <- f[, g$j[q]] - a[q]*E[,q] }
  as.vector(f)
}
relax_to_min <- function(prot, frustrated = TRUE, tol = 1e-11, maxit = 200) {
  p <- prot
  for (it in seq_len(maxit)) {
    f <- enm_force(p); if (sqrt(sum(f^2)) < tol) break
    s <- spec(hess(p, frustrated))
    Cm <- s$vec %*% ((1/s$val) * t(s$vec))
    p$nodes$xyz <- matrix(as.vector(get_xyz(p)) + as.vector(Cm %*% f), nrow = 3)
  }
  p$graph$dij <- dij_edge(matrix(as.vector(get_xyz(p)), nrow=3), get_graph(p)$i, get_graph(p)$j)
  attr(p, "fres") <- sqrt(sum(enm_force(p)^2)); attr(p, "iter") <- it
  p
}
