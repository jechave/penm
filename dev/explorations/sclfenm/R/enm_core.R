## Minimal self-contained ENM. No penm dependency: this is a testbed.

## ---- geometry ----
dists <- function(r, pr) { R <- matrix(r, nrow = 3)
  sqrt(colSums((R[, pr[,2], drop=FALSE] - R[, pr[,1], drop=FALSE])^2)) }

evecs <- function(r, pr) { R <- matrix(r, nrow = 3)
  D <- R[, pr[,2], drop=FALSE] - R[, pr[,1], drop=FALSE]
  sweep(D, 2, sqrt(colSums(D^2)), "/") }          # 3 x nedge

## ---- potential ----
## V = sum_ij [ v0_ij + 1/2 k_ij (d_ij - l_ij)^2 ]
venm <- function(r, pr, k, l, v0 = 0) sum(v0) + 0.5 * sum(k * (dists(r, pr) - l)^2)

## ---- force: F_i = sum_j k_ij (d_ij - l_ij) e_ij  (on i, toward j) ----
force <- function(r, pr, k, l) {
  N <- length(r)/3; E <- evecs(r, pr); d <- dists(r, pr)
  a <- k * (d - l); f <- matrix(0, 3, N)
  for (e in seq_len(nrow(pr))) {
    f[, pr[e,1]] <- f[, pr[e,1]] + a[e] * E[, e]
    f[, pr[e,2]] <- f[, pr[e,2]] - a[e] * E[, e]
  }
  as.vector(f)
}

## ---- Hessian ----
## frustrated = TRUE keeps the transverse term with g = l/d - 1
## (verified numerically: this sign, with bracket (ee^T - I))
hessian <- function(r, pr, k, l, frustrated = FALSE) {
  N <- length(r)/3; E <- evecs(r, pr); d <- dists(r, pr)
  K <- matrix(0, 3*N, 3*N); I3 <- diag(3)
  for (e in seq_len(nrow(pr))) {
    i <- pr[e,1]; j <- pr[e,2]
    ee <- tcrossprod(E[, e]); g <- if (frustrated) l[e]/d[e] - 1 else 0
    Kij <- -k[e] * (ee + g * (ee - I3))
    ii <- (3*i-2):(3*i); jj <- (3*j-2):(3*j)
    K[ii, jj] <- Kij; K[jj, ii] <- t(Kij)
  }
  for (i in seq_len(N)) { ii <- (3*i-2):(3*i)
    blk <- matrix(0,3,3)
    for (j in seq_len(N)) if (j != i) blk <- blk - K[ii, (3*j-2):(3*j)]
    K[ii, ii] <- blk }
  K
}

## ---- pseudo-inverse over non-rigid modes ----
nma <- function(K, tol = 1e-8) {
  e <- eigen(K, symmetric = TRUE)
  keep <- e$values > tol * max(abs(e$values))
  list(evalue = e$values[keep], umat = e$vectors[, keep, drop = FALSE],
       nzero = sum(!keep), min_eval = min(e$values))
}
cmat_of <- function(nm) nm$umat %*% ((1/nm$evalue) * t(nm$umat))

## ---- build a relaxed network on a given structure (the ENM "fit") ----
## contact set from a cutoff on the reference structure
build_relaxed <- function(r, d_max = 10.5, k_val = 1) {
  N <- length(r)/3
  pr <- t(combn(N, 2)); d <- dists(r, pr)
  keep <- d <= d_max
  pr <- pr[keep, , drop = FALSE]
  list(r = r, pr = pr, k = rep(k_val, nrow(pr)), l = dists(r, pr))  # l = d  => relaxed
}

## ---------------------------------------------------------------------------
## THE TEST PROTEIN: 2acy chain A (acylphosphatase, 98 residues), the same
## structure every penm fixture is built from. Ca coordinates only.
## ---------------------------------------------------------------------------
prot_2acy <- function() {
  e <- new.env()
  ## locate the package root by walking up from here, so this keeps working
  ## wherever the folder is moved to
  d <- normalizePath(".")
  while (!file.exists(file.path(d, "DESCRIPTION")) && dirname(d) != d) d <- dirname(d)
  f <- file.path(d, "data", "pdb_2acy_A.rda")
  stopifnot("cannot locate penm's data/pdb_2acy_A.rda" = file.exists(f))
  load(f, envir = e)
  pdb <- get("pdb_2acy_A", envir = e)
  ca  <- pdb$atom[pdb$calpha, ]
  list(r = as.vector(t(as.matrix(ca[, c("x","y","z")]))),
       resno = ca$resno, resid = ca$resid, N = nrow(ca))
}

