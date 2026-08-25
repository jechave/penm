## Phase 1 verification. Items 1-6 of the plan.
## Each assertion must be able to FAIL -- test_can_fail.R proves it does.
##
## Run: Rscript checks/test_core.R

suppressMessages(library(here))
suppressMessages(devtools::load_all(here::here("..", "..", ".."), quiet = TRUE))
source(here("R", "klfenm_core.R"))

ok <- function(label, pass, detail = "") {
  cat(sprintf("%-4s %-52s %s\n", if (pass) "OK" else "FAIL", label, detail))
  if (!pass) assign("FAILED", TRUE, envir = .GlobalEnv)
}
FAILED <- FALSE

data(pdb_2acy_A, package = "penm")
wt_prot <- set_enm(pdb_2acy_A, node = "ca", model = "anm",
                   d_max = 10.5, frustrated = FALSE)
wt <- state_from_prot(wt_prot, d_max = 10.5)
N  <- ncol(wt$R)

cat("\n== 1. Wild-type identity: the all-pairs state IS the standard ANM ==\n")

pg <- penm:::get_graph(wt_prot)   # internal: read-only use, package untouched
key_pkg <- paste(pg$i, pg$j, sep = "-")
a <- which(wt$k > 0)
key_klf <- paste(wt$pr$i[a], wt$pr$j[a], sep = "-")
ok("active set == penm's contact set",
   setequal(key_pkg, key_klf),
   sprintf("%d vs %d edges", length(key_klf), length(key_pkg)))

Kp <- matrix(as.vector(penm::get_kmat(wt_prot)), nrow = 3 * N)
Kk <- klf_hessian(wt, frustrated = TRUE)   # at wt, l == d so g == 0 anyway
ok("Hessian == penm's kmat", max(abs(Kk - Kp)) < 1e-9,
   sprintf("max|dK| = %.2e (scale %.1f)", max(abs(Kk - Kp)), max(abs(Kp))))

ok("wild type carries no strain", klf_strain(wt)$max_abs < 1e-12,
   sprintf("max|d-l| = %.2e", klf_strain(wt)$max_abs))
ok("wild type is stationary", sqrt(sum(klf_force(wt)^2)) < 1e-12,
   sprintf("|F| = %.2e", sqrt(sum(klf_force(wt)^2))))
ok("V(wt) == 0", abs(klf_v(wt)) < 1e-20, sprintf("V = %.2e", klf_v(wt)))

cat("\n== 2. The minimiser reaches a TRUE minimum ==\n")

set.seed(101)
mu <- klf_mutate(wt, site = 40, sigma = 0.3, k_update = TRUE)
fres <- attr(mu, "fres")
ok("residual force < 1e-10", fres < 1e-10,
   sprintf("|F| = %.2e in %d iters", fres, attr(mu, "iter")))

sp <- klf_spectrum(klf_hessian(mu, frustrated = TRUE))
ok("exactly 6 zero modes", sp$n_zero == 6, sprintf("n_zero = %d", sp$n_zero))
ok("no negative eigenvalues (a minimum, not a saddle)",
   min(sp$raw_value) > -1e-8,
   sprintf("lowest raw eigenvalue = %+.3e", min(sp$raw_value)))
ok("mutant IS strained (else nothing is being tested)",
   klf_strain(mu)$max_abs > 1e-3,
   sprintf("max|d-l| = %.4f", klf_strain(mu)$max_abs))

cat("\n== 3. Exact V vs the linear-response two-term formula ==\n")
## They must agree as sigma -> 0 and DIVERGE as sigma grows. If they agree
## everywhere, the minimiser is not doing anything the LRA could not.
##
## k_update = FALSE here ON PURPOSE. With k(l) live, a perturbation that pushes
## l across d_max deletes a spring and its stored strain, so exact V and the
## two-term formula differ by a whole contact -- a real effect of the model, but
## NOT the truncation error this check is about. Measured: at sigma = 0.02 that
## already happens (a pair at l = 10.4999 -> 10.5113), giving a 26% gap that has
## nothing to do with linear response. Freezing k isolates the truncation.
C_wt <- klf_cmat(klf_hessian(wt, frustrated = FALSE))
two_term <- function(st0, dl_idx, dl) {
  ## V_stress - V_relax, with f = -k dl along e, using the WT compliance
  pr <- st0$pr
  R  <- st0$R
  D  <- R[, pr$j[dl_idx], drop = FALSE] - R[, pr$i[dl_idx], drop = FALSE]
  d  <- sqrt(colSums(D^2)); E <- sweep(D, 2, d, "/")
  f  <- matrix(0, 3, ncol(R))
  fij <- -st0$k[dl_idx] * dl
  for (q in seq_along(dl_idx)) {
    f[, pr$i[dl_idx][q]] <- f[, pr$i[dl_idx][q]] + fij[q] * E[, q]
    f[, pr$j[dl_idx][q]] <- f[, pr$j[dl_idx][q]] - fij[q] * E[, q]
  }
  fv <- as.vector(f)
  dr <- as.vector(C_wt %*% fv)
  0.5 * sum(st0$k[dl_idx] * dl^2) - 0.5 * sum(dr * as.vector(klf_hessian(st0, frustrated = FALSE) %*% dr))
}
sel40 <- which((wt$pr$i == 40 | wt$pr$j == 40) & wt$k > 0)
rel <- sapply(c(0.02, 0.1, 0.4), function(sg) {
  set.seed(7)
  dl <- rnorm(length(sel40), 0, sg)
  s2 <- wt; s2$l[sel40] <- s2$l[sel40] + dl
  ## k frozen: no contact may break, so the only gap is the LRA truncation
  s2 <- klf_minimise(s2)
  100 * abs(klf_v(s2) - two_term(wt, sel40, dl)) / klf_v(s2)
})
cat(sprintf("     sigma 0.02 / 0.10 / 0.40  ->  gap %.3f%% / %.3f%% / %.3f%%\n", rel[1], rel[2], rel[3]))
ok("exact and LRA agree at small sigma", rel[1] < 0.5, sprintf("%.3f%%", rel[1]))
ok("and diverge at large sigma", rel[3] > 3 * rel[1],
   sprintf("%.2f%% vs %.3f%%", rel[3], rel[1]))

cat("\n== 4. k = k(l) consistency after mutation ==\n")
kf <- k_fun_of(mu$k_model)
k_expect <- do.call(kf, c(list(lij = mu$l, sdij = mu$sdij), mu$k_par))
ok("k_ij == k(l_ij) for every pair", identical(as.numeric(mu$k), as.numeric(k_expect)),
   sprintf("max|dk| = %.2e", max(abs(mu$k - k_expect))))

## and with k_update = FALSE, k must NOT follow l
set.seed(101)
mu0 <- klf_mutate(wt, site = 40, sigma = 0.3, k_update = FALSE)
k_would <- do.call(kf, c(list(lij = mu0$l, sdij = mu0$sdij), mu0$k_par))
ok("k_update=FALSE freezes k (active set unchanged)",
   n_active(mu0) == n_active(wt),
   sprintf("%d active, vs %d if k had followed l", n_active(mu0), sum(k_would > 0)))

cat("\n== 5. Contact events occur, in BOTH directions ==\n")
## Both directions matters: if only breaking happens, the network monotonically
## thins and the k(l) mechanism is half-dead. That is what an earlier version
## did, by perturbing only a site's ACTIVE pairs -- 14 broken, 0 formed.
set.seed(2024)
st <- wt; nb <- nf <- 0L
for (s in 1:25) {
  st <- klf_mutate(st, sample(N, 1), sigma = 0.3, k_update = TRUE, radius = 12)
  nb <- nb + attr(st, "n_broken"); nf <- nf + attr(st, "n_formed")
}
ok("contacts break along a walk", nb > 0, sprintf("%d broken", nb))
ok("contacts also FORM along a walk", nf > 0, sprintf("%d formed", nf))
cat(sprintf("     active %d -> %d over 25 mutations\n", n_active(wt), n_active(st)))

## and the perturbation really does reach inactive pairs
sel_all <- which((wt$pr$i == 40 | wt$pr$j == 40) & wt$sdij >= 1 &
                 (pair_dist(wt$R, wt$pr) <= 12 | wt$k > 0))
ok("mutation reaches inactive pairs", sum(wt$k[sel_all] == 0) > 0,
   sprintf("%d of %d perturbed pairs at site 40 are inactive",
           sum(wt$k[sel_all] == 0), length(sel_all)))

cat("\n== 6. Reversibility of the parameters ==\n")
set.seed(55)
sel <- which((wt$pr$i == 11 | wt$pr$j == 11) & wt$k > 0)
dl  <- rnorm(length(sel), 0, 0.3)
f1 <- wt; f1$l[sel] <- f1$l[sel] + dl; f1 <- refresh_k(f1); f1 <- klf_minimise(f1)
b1 <- f1; b1$l[sel] <- b1$l[sel] - dl; b1 <- refresh_k(b1); b1 <- klf_minimise(b1)
ok("l returns exactly", max(abs(b1$l - wt$l)) < 1e-12,
   sprintf("max|dl| = %.2e", max(abs(b1$l - wt$l))))
ok("k returns exactly", identical(as.numeric(b1$k), as.numeric(wt$k)))
## Compare INTERNAL geometry: the Newton step uses the pseudo-inverse, which is
## orthogonal to the six rigid-body modes but does not pin the frame, so the
## molecule drifts bodily while its internal structure returns exactly. Measured
## raw max|dR| = 0.41 vs 1.2e-14 after superposition -- all of it rigid-body.
kabsch_rmsd <- function(X, Y) {
  A <- X - rowMeans(X); B <- Y - rowMeans(Y)
  s <- svd(A %*% t(B)); dsign <- sign(det(s$v %*% t(s$u)))
  Brot <- t(s$v %*% diag(c(1, 1, dsign)) %*% t(s$u)) %*% B
  list(rmsd = sqrt(mean(colSums((Brot - A)^2))), max = max(abs(Brot - A)))
}
kb <- kabsch_rmsd(wt$R, b1$R)
ok("internal structure returns (after superposition)", kb$max < 1e-9,
   sprintf("max|dR| = %.2e, RMSD = %.2e (raw, unsuperposed: %.3f)",
           kb$max, kb$rmsd, max(abs(b1$R - wt$R))))
ok("and V returns to zero", abs(klf_v(b1)) < 1e-20,
   sprintf("V = %.2e", klf_v(b1)))

cat("\n")
if (FAILED) { cat("SOME CHECKS FAILED\n"); quit(status = 1) } else cat("all checks passed\n")
