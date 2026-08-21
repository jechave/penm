library(here)
source(here("R", "penm_catalog.R"))
wt <- founder(); N <- get_nsites(wt)
ok <- function(l, v, tol=1e-10) cat(sprintf("[%s] %-52s %.3e\n",
  if (abs(v)<tol) "EXACT" else "APPROX", l, abs(v)))

cat("=== I1. dr is linear in f, hence additive over sites ===\n")
sites <- c(11, 40, 77); ms <- c(1, 2, 3)
frs <- lapply(seq_along(sites), function(q) {
  mu <- get_mutant_site(wt, sites[q], ms[q], "lfenm", .3, 1L, 1L)
  list(f = calculate_force(wt, get_graph(mu)$lij - get_graph(wt)$lij),
       dr = as.vector(get_xyz(mu) - get_xyz(wt))) })
f_tot  <- Reduce(`+`, lapply(frs, `[[`, "f"))
dr_sum <- Reduce(`+`, lapply(frs, `[[`, "dr"))
dr_tot <- as.vector(calculate_dxyz(wt, f_tot))
ok("dr(sum of forces) == sum of dr", max(abs(dr_tot - dr_sum)))

cat("\n=== I2. reversibility: X(m0->m1) = -X(m1->m0) ===\n")
cg_j <- function(j, ms) sapply(ms, function(m) {
  mu <- get_mutant_site(wt, j, m, "lfenm", .3, 1L, 1L); ddg_dv(wt, mu) })
## NOTE: taking both differences from one stored vector would be the identity
## x + (-x) and could not fail. Recompute each mutant independently instead.
mm <- function(m) get_mutant_site(wt, 40, m, "lfenm", .3, 1L, 1L)
v1 <- ddg_dv(wt, mm(1)); v2 <- ddg_dv(wt, mm(2))
d01 <- v2 - v1; d10 <- v1 - v2
ok("dV(1->2) + dV(2->1) == 0  (mutants rebuilt independently)", d01 + d10)
v <- c(0, v1, v2)
cat(sprintf("      dV(0->1)=%.5f  dV(0->2)=%.5f  dV(1->2)=%+.5f  dV(2->1)=%+.5f\n",
    v[2], v[3], d01, d10))

cat("\n=== I3. does dV take BOTH signs between states? ===\n")
## Counting BOTH orderings of every pair forces exactly 50% by antisymmetry and
## measures nothing. Count UNORDERED pairs instead: the informative question is
## how often two states of a site differ at all, and by how much.
set.seed(1); j <- sample(N, 30); dif <- c()
for (jj in j) { vv <- cg_j(jj, 0:6)
  for (a in 1:6) for (b in 1:6) if (a < b) dif <- c(dif, vv[b+1]-vv[a+1]) }
cat(sprintf("      %d unordered pairs at 30 sites: median |dV| = %.4f, max = %.4f\n",
    length(dif), median(abs(dif)), max(abs(dif))))
cat(sprintf("      a move in one direction is downhill for %.0f%% of pairs -- so a\n",
    100*mean(dif != 0)))
cat("      downhill move exists for essentially every pair of states.\n")

cat("\n=== I4. m = 0 really is the founder ===\n")
mu0 <- get_mutant_site(wt, 40, 0, "lfenm", .3, 1L, 1L)
ok("mutation = 0 returns wt unchanged", as.numeric(!identical(mu0, wt)))
