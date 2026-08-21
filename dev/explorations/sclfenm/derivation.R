## ============================================================================
## Every step of the derivation, checked numerically on 2acy chain A.
## NOTATION (fixed, used everywhere):
##
##   A Hamiltonian is fixed by its spring parameters. Three appear:
##     H_wt    : l = l_wt                    (wild type)
##     H_mut   : l = l_wt + dl               (LFENM mutant; KEEPS the strain)
##     H_ref   : l = d^e(r^e_mut)            (network REFIT on the mutant structure)
##
##   V_X(r)   value of Hamiltonian X at conformation r          [a function]
##   r^e_X    equilibrium conformation of X, where F_X(r^e_X)=0
##   V^min_X  = V_X(r^e_X)                                       [a number]
##
##   d_ij(r)  distance in conformation r
##   d^e_ij   = d_ij(r^e) for the relevant Hamiltonian
## ============================================================================
source("sclfenm.R"); source("fig_setup.R")

p  <- prot_2acy()
wt <- new_state(p$r, d_max = 10.5)
K_wt <- state_K(wt); nm_wt <- nma(K_wt); C_wt <- cmat_of(nm_wt)
cat(sprintf("2acy chain A: N=%d, edges=%d, zero modes=%d\n\n", p$N, nrow(wt$pr), nm_wt$nzero))

ok <- function(lbl, lhs, rhs, tol = 1e-8) {
  d <- abs(lhs - rhs); rel <- d/max(abs(rhs), 1e-12)
  cat(sprintf("[%s] %-58s %12.6f vs %12.6f  (err %.2e)\n",
      if (d < tol || rel < tol) "OK  " else "DIFF", lbl, lhs, rhs, d))
  invisible(d)
}

## --- the mutation -----------------------------------------------------------
site <- 40
mask <- wt$pr[,1]==site | wt$pr[,2]==site
set.seed(11); dl <- rep(0, nrow(wt$pr)); dl[mask] <- rnorm(sum(mask), 0, 0.3)
l_wt  <- wt$l
l_mut <- l_wt + dl
cat(sprintf("mutating site %d (%d contacts), sigma=0.3\n\n", site, sum(mask)))

## --- STEP 0. the wild type is relaxed: l_wt = d^e(r^e_wt) -------------------
cat("STEP 0 -- the wild type\n")
ok("l_wt == d(r^e_wt)  (max abs deviation)",
   max(abs(dists(wt$r, wt$pr) - l_wt)), 0)
ok("V^min_wt = V_wt(r^e_wt) = 0", venm(wt$r, wt$pr, wt$k, l_wt), 0)
ok("|F_wt(r^e_wt)| = 0", sqrt(sum(force(wt$r, wt$pr, wt$k, l_wt)^2)), 0)
cat("\n")

## --- STEP 1. H_mut evaluated at r^e_wt: the STRESS energy -------------------
## V_mut(r^e_wt) = 1/2 sum k (d^e_wt - l_wt - dl)^2 = 1/2 sum k dl^2
cat("STEP 1 -- H_mut at the WILD-TYPE conformation (stress)\n")
V_mut_at_wt <- venm(wt$r, wt$pr, wt$k, l_mut)
ok("V_mut(r^e_wt) == 1/2 sum k dl^2", V_mut_at_wt, 0.5*sum(wt$k*dl^2))
cat("\n")

## --- STEP 2. the induced force and the linear response ----------------------
## F_mut(r^e_wt)_i = sum_j k(d-l_mut) e_ij, and with d = l_wt this is -k dl
cat("STEP 2 -- force induced at r^e_wt, and the linear response\n")
f  <- force(wt$r, wt$pr, wt$k, l_mut)
ok("net force sums to zero (translational invariance)", sum(abs(rowSums(matrix(f,3)))), 0, 1e-9)
## torque about the centroid
R <- matrix(wt$r,3); Fm <- matrix(f,3); cen <- rowMeans(R)
tau <- rowSums(sapply(seq_len(p$N), function(i) {
  v <- R[,i]-cen; c(v[2]*Fm[3,i]-v[3]*Fm[2,i], v[3]*Fm[1,i]-v[1]*Fm[3,i],
                    v[1]*Fm[2,i]-v[2]*Fm[1,i]) }))
ok("net torque is zero (central forces) -> f is in range(K)", sum(abs(tau)), 0, 1e-8)
dr <- as.vector(C_wt %*% f)
ok("K dr == f  (dr = C f solves the linear system)", max(abs(K_wt%*%dr - f)), 0, 1e-8)
r_mut_lra <- wt$r + dr
cat("\n")

## --- STEP 3. is r^e_mut = r^e_wt + dr ? Only to O(dl^2). --------------------
cat("STEP 3 -- the mutant's equilibrium conformation\n")
st_mut <- list(r = r_mut_lra, pr = wt$pr, k = wt$k, l = l_mut, frustrated = FALSE)
rq <- state_req(st_mut)
cat(sprintf("     |F_mut| at r^e_wt + dr (LRA)     = %.3e\n",
    sqrt(sum(force(r_mut_lra, wt$pr, wt$k, l_mut)^2))))
cat(sprintf("     |F_mut| at true r^e_mut          = %.3e   (%d iterations)\n",
    rq$fres, rq$iter))
cat(sprintf("     |r_LRA - r^e_mut|                = %.3e\n", sqrt(sum((r_mut_lra-rq$r)^2))))
r_mut <- rq$r      # use the TRUE minimum from here on
cat("\n")

## --- STEP 4. V^min_mut: the mutant Hamiltonian at ITS OWN minimum -----------
cat("STEP 4 -- V^min_mut (this is what dV compares; it is NOT zero)\n")
V_min_mut <- venm(r_mut, wt$pr, wt$k, l_mut)
dV_exact  <- V_min_mut - venm(wt$r, wt$pr, wt$k, l_wt)
cat(sprintf("     V^min_mut = %.6f\n", V_min_mut))
cat(sprintf("     dV = V^min_mut - V^min_wt = %.6f\n", dV_exact))
## the two-term identity, with dr measured to the TRUE minimum
dr_true <- r_mut - wt$r
V_stress <- 0.5*sum(wt$k*dl^2)
V_relax  <- 0.5*as.numeric(crossprod(dr_true, K_wt %*% dr_true))
ok("dV == V_stress - V_relax   (to O(dl^3))", V_stress - V_relax, dV_exact, 1e-3)
cat(sprintf("     V_stress = %.6f   V_relax = %.6f   err = %.2e\n\n",
    V_stress, V_relax, (V_stress-V_relax)-dV_exact))

## --- STEP 5. H_ref: the REFIT Hamiltonian. This is the one with V^min = 0 ---
cat("STEP 5 -- H_ref, the refit network (a DIFFERENT Hamiltonian)\n")
ref <- build_relaxed(r_mut, wt$d_max, wt$k_val)
ok("l_ref == d(r^e_mut)  (refit springs are relaxed)",
   max(abs(dists(ref$r, ref$pr) - ref$l)), 0)
ok("V^min_ref = V_ref(r^e_mut) = 0", venm(ref$r, ref$pr, ref$k, ref$l), 0)
cat(sprintf("     edges: wt %d -> refit %d\n", nrow(wt$pr), nrow(ref$pr)))
cat(sprintf("     BUT V^min_mut = %.6f. The two Hamiltonians share the structure\n", V_min_mut))
cat( "     r^e_mut and disagree about its energy. That difference is what must\n")
cat( "     be carried in v_off.\n\n")

## --- STEP 6. the offset makes the refit energy agree with H_mut ------------
cat("STEP 6 -- with the offset\n")
ok("V^min_ref + v_off == V^min_mut", venm(ref$r,ref$pr,ref$k,ref$l) + dV_exact, V_min_mut)
