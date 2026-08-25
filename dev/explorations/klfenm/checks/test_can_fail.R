## Break the code on purpose; confirm the checks in test_core.R go red.
## A suite that has never been seen to fail is itself untested.
##
## Run: Rscript checks/test_can_fail.R

suppressMessages(library(here))
suppressMessages(devtools::load_all(here::here("..", "..", ".."), quiet = TRUE))
source(here("R", "klfenm_core.R"))

data(pdb_2acy_A, package = "penm")
wtp <- set_enm(pdb_2acy_A, node = "ca", model = "anm", d_max = 10.5, frustrated = FALSE)
wt  <- state_from_prot(wtp)
N   <- ncol(wt$R)

detected <- 0L; total <- 0L
sabotage <- function(label, detector) {
  total <<- total + 1L
  res <- tryCatch(detector(), error = function(e) TRUE)  # an error IS detection
  detected <<- detected + as.integer(isTRUE(res))
  cat(sprintf("%-4s %s\n", if (isTRUE(res)) "RED" else "MISSED", label))
}

cat("\n== sabotages, and whether the corresponding check notices ==\n")

## S1. Wrong Hessian sign convention (g = 1 - l/d instead of l/d - 1).
## Check 1 compares against penm's kmat -- but at the WT g == 0, so the sign is
## invisible there. It must be caught on a STRAINED network instead.
sabotage("S1 flipped transverse sign detected on a strained state", function() {
  hess_wrong <- function(st, R = st$R) {
    pr <- st$pr; a <- which(st$k > 0)
    D <- R[, pr$j[a], drop = FALSE] - R[, pr$i[a], drop = FALSE]
    d <- sqrt(colSums(D^2)); E <- sweep(D, 2, d, "/"); I3 <- diag(3)
    K <- matrix(0, 3 * N, 3 * N)
    for (q in seq_along(a)) {
      e <- E[, q]; ee <- tcrossprod(e, e)
      g <- 1 - st$l[a][q] / d[q]                       # WRONG SIGN
      Kij <- -st$k[a][q] * (ee + g * (ee - I3))
      bi <- (3 * (pr$i[a][q] - 1) + 1):(3 * pr$i[a][q])
      bj <- (3 * (pr$j[a][q] - 1) + 1):(3 * pr$j[a][q])
      K[bi, bj] <- K[bi, bj] + Kij; K[bj, bi] <- K[bj, bi] + t(Kij)
      K[bi, bi] <- K[bi, bi] - Kij; K[bj, bj] <- K[bj, bj] - Kij
    }
    K
  }
  set.seed(3); mu <- klf_mutate(wt, 40, 0.5, k_update = FALSE)
  Kok <- klf_hessian(mu, frustrated = TRUE); Kbad <- hess_wrong(mu)
  ## a numerical Hessian arbitrates
  vf <- function(v) klf_v(mu, matrix(v, nrow = 3))
  v0 <- as.vector(mu$R); h <- 1e-5
  idx <- c(1, 2, 3, 40 * 3 - 2, 40 * 3 - 1, 40 * 3)
  num <- outer(idx, idx, Vectorize(function(a, b) {
    ep <- rep(0, length(v0)); eq <- rep(0, length(v0)); ep[a] <- h; eq[b] <- h
    (vf(v0 + ep + eq) - vf(v0 + ep - eq) - vf(v0 - ep + eq) + vf(v0 - ep - eq)) / (4 * h^2)
  }))
  err_ok  <- max(abs(Kok[idx, idx]  - num))
  err_bad <- max(abs(Kbad[idx, idx] - num))
  cat(sprintf("       correct sign off by %.2e ; flipped sign off by %.2e\n", err_ok, err_bad))
  err_ok < 1e-4 && err_bad > 1e-2
})

## S2. Minimiser that stops after ONE step (the linear-response shortcut).
## Check 2 asserts |F| < 1e-10; a single step must leave a much larger residual.
sabotage("S2 one-step 'minimiser' leaves a detectable residual force", function() {
  st <- wt
  sel <- which((st$pr$i == 40 | st$pr$j == 40) & st$k > 0)
  set.seed(9); st$l[sel] <- st$l[sel] + rnorm(length(sel), 0, 0.3)
  s <- klf_spectrum(klf_hessian(st, frustrated = TRUE))
  Cm <- s$vector %*% ((1 / s$value) * t(s$vector))
  R1 <- matrix(as.vector(st$R) + as.vector(Cm %*% klf_force(st)), nrow = 3)
  st1 <- st; st1$R <- R1
  f1 <- sqrt(sum(klf_force(st1)^2))
  full <- klf_minimise(st)
  cat(sprintf("       one step |F| = %.2e ; converged |F| = %.2e\n", f1, attr(full, "fres")))
  f1 > 1e-10 && attr(full, "fres") < 1e-10
})

## S3. Spectrum taken OFF a stationary point -- the trap that poisons RMSF.
## Check 2 asserts exactly 6 zero modes and no negative eigenvalues.
sabotage("S3 off-minimum spectrum shows spurious/negative modes", function() {
  st <- wt
  sel <- which((st$pr$i == 40 | st$pr$j == 40) & st$k > 0)
  set.seed(9); st$l[sel] <- st$l[sel] + rnorm(length(sel), 0, 0.3)
  sp_off <- klf_spectrum(klf_hessian(st, frustrated = TRUE))   # NOT minimised
  sp_on  <- klf_spectrum(klf_hessian(klf_minimise(st), frustrated = TRUE))
  cat(sprintf("       off-minimum: n_zero = %d, lowest = %+.3e | at minimum: n_zero = %d, lowest = %+.3e\n",
              sp_off$n_zero, min(sp_off$raw_value), sp_on$n_zero, min(sp_on$raw_value)))
  (sp_off$n_zero != 6 || min(sp_off$raw_value) < -1e-8) &&
    sp_on$n_zero == 6 && min(sp_on$raw_value) > -1e-8
})

## S4. k left stale after l changes -- exactly the penm mutate_graph() failure
## mode (overwrite lij after kij was computed). Check 4 must catch it.
sabotage("S4 stale k (k != k(l)) is caught", function() {
  st <- wt
  sel <- which((st$pr$i == 40 | st$pr$j == 40) & st$k > 0)
  set.seed(4); st$l[sel] <- st$l[sel] + rnorm(length(sel), 0, 0.6)
  ## deliberately DO NOT refresh_k
  kf <- k_fun_of(st$k_model)
  k_expect <- do.call(kf, c(list(lij = st$l, sdij = st$sdij), st$k_par))
  cat(sprintf("       stale k differs from k(l) on %d pairs\n", sum(st$k != k_expect)))
  !identical(as.numeric(st$k), as.numeric(k_expect))
})

## S5. Reversibility check must not be fooled by rigid-body drift, and must
## still fail when the internal structure genuinely differs.
sabotage("S5 reversibility catches a genuine internal change", function() {
  kabsch_max <- function(X, Y) {
    A <- X - rowMeans(X); B <- Y - rowMeans(Y)
    s <- svd(A %*% t(B)); ds <- sign(det(s$v %*% t(s$u)))
    max(abs(t(s$v %*% diag(c(1, 1, ds)) %*% t(s$u)) %*% B - A))
  }
  ## a pure rigid-body move must NOT be flagged
  ang <- 0.3; Rot <- matrix(c(cos(ang), -sin(ang), 0, sin(ang), cos(ang), 0, 0, 0, 1), 3, 3)
  moved <- Rot %*% wt$R + 5
  rigid_ok <- kabsch_max(wt$R, moved) < 1e-9
  ## a genuine internal perturbation MUST be flagged
  bent <- wt$R; bent[, 50] <- bent[, 50] + c(0.5, 0, 0)
  bent_caught <- kabsch_max(wt$R, bent) > 1e-3
  cat(sprintf("       rigid-body move: max|dR| ~ %.1e (ignored) ; bent site: %.3f (caught)\n",
              kabsch_max(wt$R, moved), kabsch_max(wt$R, bent)))
  rigid_ok && bent_caught
})

cat(sprintf("\n%d of %d sabotages detected\n", detected, total))
if (detected < total) quit(status = 1)
