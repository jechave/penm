# Run from the package root: Rscript dev/explorations/motion_comparisons/complement.R
#
# Does a network share its eigenvectors with its complement network?
# GNM (graph Laplacian): yes, exactly. ANM (3D Hessian): no.
suppressMessages(devtools::load_all(quiet = TRUE))
wt <- set_enm(pdb_2acy_A, node = "ca", model = "anm", d_max = 10.5)
xyz <- matrix(get_xyz(wt), nrow = 3)
n <- ncol(xyz)
d <- unname(as.matrix(dist(t(xyz))))

# GNM: Kirchhoff matrix of the contact graph and of its complement
a <- (d <= 10.5) * 1; diag(a) <- 0
lap <- function(a) diag(rowSums(a)) - a
g <- lap(a)
g_bar <- lap(1 - a - diag(n))
cat("GNM: L + L_bar == n I - J:", isTRUE(all.equal(g + g_bar, n * diag(n) - matrix(1, n, n))), "\n")
cat("GNM: commutator norm ||L L_bar - L_bar L|| =", norm(g %*% g_bar - g_bar %*% g, "F"), "\n")
lg <- eigen(g, symmetric = TRUE); lgb <- eigen(g_bar, symmetric = TRUE)
cat("GNM: eigenvalues of L_bar == n - eigenvalues of L (non-zero modes, reversed):",
    isTRUE(all.equal(sort(lgb$values[1:(n - 1)]), sort(n - lg$values[1:(n - 1)]))), "\n")

# ANM: Hessian of the contact network and of its complement (unit springs)
hess <- function(a) {
  h <- matrix(0, 3 * n, 3 * n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) if (a[i, j] != 0) {
    e <- (xyz[, j] - xyz[, i]) / d[i, j]
    b <- a[i, j] * tcrossprod(e)
    ii <- 3 * i - 2:0; jj <- 3 * j - 2:0
    h[ii, jj] <- h[ii, jj] - b; h[jj, ii] <- h[jj, ii] - b
    h[ii, ii] <- h[ii, ii] + b; h[jj, jj] <- h[jj, jj] + b
  }
  h
}
h <- hess(a)
h_bar <- hess(1 - a - diag(n))
cat("ANM: check h equals penm kmat:", isTRUE(all.equal(h, unname(as.matrix(get_kmat(wt))))), "\n")
cat("ANM: relative commutator ||H H_bar - H_bar H|| / (||H|| ||H_bar||) =",
    norm(h %*% h_bar - h_bar %*% h, "F") / (norm(h, "F") * norm(h_bar, "F")), "\n")
# overlap of each H mode with its best H_bar mode
u <- eigen(h, symmetric = TRUE)$vectors[, 1:(3 * n - 6)]
v <- eigen(h_bar, symmetric = TRUE)$vectors[, 1:(3 * n - 6)]
best <- apply(abs(crossprod(u, v)), 1, max)
cat("ANM: best |overlap| of H modes with H_bar modes, quantiles:", round(quantile(best, c(0, .25, .5, .75, 1)), 3), "\n")
