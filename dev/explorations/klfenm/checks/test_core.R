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
wt <- klfenm_set_state(wt_prot, d_max = 10.5)
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
Kk <- klfenm_kmat(wt, frustrated = TRUE)   # at wt, l == d so g == 0 anyway
ok("Hessian == penm's kmat", max(abs(Kk - Kp)) < 1e-9,
   sprintf("max|dK| = %.2e (scale %.1f)", max(abs(Kk - Kp)), max(abs(Kp))))

ok("wild type carries no strain", klfenm_strain(wt)$max_abs < 1e-12,
   sprintf("max|d-l| = %.2e", klfenm_strain(wt)$max_abs))
ok("wild type is stationary", sqrt(sum(energy_gradient(wt)^2)) < 1e-12,
   sprintf("|F| = %.2e", sqrt(sum(energy_gradient(wt)^2))))
ok("V(wt) == 0", abs(klfenm_energy(wt)) < 1e-20, sprintf("V = %.2e", klfenm_energy(wt)))

cat("\n== 2. The minimiser reaches a TRUE minimum ==\n")

set.seed(101)
mu <- klfenm_mutate_site(wt, site = 40, sigma = 0.3, k_update = TRUE)
fres <- attr(mu, "fres")
ok("residual force < 1e-10", fres < 1e-10,
   sprintf("|F| = %.2e in %d iters", fres, attr(mu, "iter")))

sp <- klfenm_nma(klfenm_kmat(mu, frustrated = TRUE))
ok("exactly 6 zero modes", sp$n_zero == 6, sprintf("n_zero = %d", sp$n_zero))
ok("no negative eigenvalues (a minimum, not a saddle)",
   min(sp$raw_value) > -1e-8,
   sprintf("lowest raw eigenvalue = %+.3e", min(sp$raw_value)))
ok("mutant IS strained (else nothing is being tested)",
   klfenm_strain(mu)$max_abs > 1e-3,
   sprintf("max|d-l| = %.4f", klfenm_strain(mu)$max_abs))

cat("\n== 3. The reference's own energy is not silently dropped ==\n")
## NOT "the split is exact": stress + relax = min telescopes --
##   [Vm(Rref) - Vref(Rref)] + [Vm(Rm) - Vm(Rref)] = Vm(Rm) - Vref(Rref)
## -- so it holds for any three numbers of that form and verifies nothing. It is
## computed below only to confirm the implementations match their definitions.
##
## The assertion with content is the last one: V_mut(r_ref), which is what the
## sloppy notation "V_stress" invites you to write, equals dV_stress ONLY at a
## relaxed reference. At a strained one they differ by exactly V_ref -- measured
## sevenfold at 15 substitutions. That can fail, and would if the reference's
## energy were being dropped.

check_split <- function(ref, label) {
  set.seed(7)
  sel <- which((ref$pr$i == 40 | ref$pr$j == 40) & ref$sdij >= 1)
  m <- ref; m$l[sel] <- m$l[sel] + rnorm(length(sel), 0, 0.3)
  m <- refresh_k(m); m <- klfenm_minimise(m)
  s <- klfenm_delta_v_stress(ref, m); r <- klfenm_delta_v_relax(ref, m); tot <- klfenm_delta_v_min(ref, m)
  cat(sprintf("     %-22s V_ref(r_ref) = %8.4f | dV_stress %+8.4f  dV_relax %+8.4f  dV %+8.4f\n",
              label, klfenm_energy(ref), s, r, tot))
  list(gap = abs((s + r) - tot), relax = r, stress = s,
       vref = klfenm_energy(ref), naive = klfenm_energy(m, R = ref$R))
}

a <- check_split(wt, "relaxed reference")
ok("implementations match their definitions (telescoping, cannot fail)",
   a$gap < 1e-12, sprintf("gap = %.2e", a$gap))
ok("dV_relax <= 0 (relaxation returns energy)", a$relax <= 1e-12,
   sprintf("%.4f", a$relax))

## now a strained reference, 15 substitutions along
set.seed(31); ref2 <- wt
for (s in 1:15) ref2 <- klfenm_mutate_site(ref2, sample(N, 1), sigma = 0.3, k_update = TRUE)
b <- check_split(ref2, "strained reference")
ok("same at a strained reference (also cannot fail)", b$gap < 1e-12,
   sprintf("gap = %.2e", b$gap))
ok("dV_relax <= 0 (strained ref)", b$relax <= 1e-12, sprintf("%.4f", b$relax))

## the trap itself: V_mut(r_ref) is NOT dV_stress once the reference is strained
ok("V_mut(r_ref) == dV_stress ONLY at a relaxed reference",
   abs(a$naive - a$stress) < 1e-12 && abs(b$naive - b$stress) > 1,
   sprintf("relaxed: identical | strained: %.4f vs %.4f, off by %.4f (= V_ref)",
           b$naive, b$stress, b$naive - b$stress))

cat("\n== 4. k = k(l) consistency after mutation ==\n")
kf <- k_fun_of(mu$k_model)
k_expect <- do.call(kf, c(list(lij = mu$l, sdij = mu$sdij), mu$k_par))
ok("k_ij == k(l_ij) for every pair", identical(as.numeric(mu$k), as.numeric(k_expect)),
   sprintf("max|dk| = %.2e", max(abs(mu$k - k_expect))))

## and with k_update = FALSE, k must NOT follow l
set.seed(101)
mu0 <- klfenm_mutate_site(wt, site = 40, sigma = 0.3, k_update = FALSE)
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
  st <- klfenm_mutate_site(st, sample(N, 1), sigma = 0.3, k_update = TRUE, radius = 12)
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
f1 <- wt; f1$l[sel] <- f1$l[sel] + dl; f1 <- refresh_k(f1); f1 <- klfenm_minimise(f1)
b1 <- f1; b1$l[sel] <- b1$l[sel] - dl; b1 <- refresh_k(b1); b1 <- klfenm_minimise(b1)
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
ok("and V returns to zero", abs(klfenm_energy(b1)) < 1e-20,
   sprintf("V = %.2e", klfenm_energy(b1)))

cat("\n")
if (FAILED) { cat("SOME CHECKS FAILED\n"); quit(status = 1) } else cat("all checks passed\n")
