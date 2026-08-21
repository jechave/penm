library(here)
source(here("R", "enm_core.R"))
set.seed(1)

## a compact random "protein": N nodes on a perturbed helix-ish blob
make_prot <- function(N = 40) {
  t <- seq(0, 6*pi, length.out = N)
  r <- rbind(4*cos(t) + rnorm(N,0,.4), 4*sin(t) + rnorm(N,0,.4), 1.2*t/2 + rnorm(N,0,.4))
  as.vector(r)
}

r_wt <- make_prot(40)
wt <- build_relaxed(r_wt, d_max = 10.5)
cat("N =", length(r_wt)/3, " edges =", nrow(wt$pr), "\n")

K_wt <- hessian(wt$r, wt$pr, wt$k, wt$l, frustrated = FALSE)
nm_wt <- nma(K_wt); C_wt <- cmat_of(nm_wt)
cat("zero modes:", nm_wt$nzero, " (expect 6)\n\n")

## ---- mutate: perturb l on edges of one site ----
site <- 11
mask <- wt$pr[,1] == site | wt$pr[,2] == site
dl <- rep(0, nrow(wt$pr)); dl[mask] <- rnorm(sum(mask), 0, 0.3)
cat("mutated site", site, "with", sum(mask), "edges\n\n")

l_mut <- wt$l + dl

## Force induced at the WT structure by the perturbation.
## F(r_wt) with l_mut: since d(r_wt) = l_wt, we get a = k*(l_wt - l_mut) = -k*dl
f <- force(r_wt, wt$pr, wt$k, l_mut)
dr <- C_wt %*% f
r_mut <- r_wt + as.vector(dr)

## ---- THE IDENTITY ----
## (A) direct: energy of the perturbed (frustrated) Hamiltonian at its relaxed structure
V_direct <- venm(r_mut, wt$pr, wt$k, l_mut)

## (B) decomposition: stress at wt structure  minus  relaxation recovered
V_stress <- venm(r_wt, wt$pr, wt$k, l_mut)          # = 1/2 sum k dl^2
V_relax  <- 0.5 * as.numeric(t(dr) %*% K_wt %*% dr)
V_decomp <- V_stress - V_relax

cat("V_stress  (1/2 sum k dl^2)     =", V_stress, "\n")
cat("check 1/2 sum k dl^2           =", 0.5*sum(wt$k*dl^2), "\n")
cat("V_relax   (1/2 dr' K dr)       =", V_relax, "\n")
cat("V_decomp  = stress - relax     =", V_decomp, "\n")
cat("V_direct  (exact, at r_mut)    =", V_direct, "\n")
cat("difference                     =", V_direct - V_decomp, "\n")
cat("relative                       =", (V_direct - V_decomp)/V_direct, "\n\n")

## Is r_mut actually the minimum of the perturbed Hamiltonian?
cat("|F| at r_wt  (perturbed H) =", sqrt(sum(force(r_wt, wt$pr, wt$k, l_mut)^2)), "\n")
cat("|F| at r_mut (perturbed H) =", sqrt(sum(force(r_mut, wt$pr, wt$k, l_mut)^2)), "\n")
