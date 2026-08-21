source("sclfenm.R"); source("fig_setup.R")
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5)
K <- state_K(wt); C <- cmat_of(nma(K))
site <- 11; mask <- wt$pr[,1]==site | wt$pr[,2]==site
set.seed(11); dl <- rep(0,nrow(wt$pr)); dl[mask] <- rnorm(sum(mask),0,.3)
l_mut <- wt$l + dl
f  <- force(wt$r, wt$pr, wt$k, l_mut)
dr <- as.vector(C %*% f); r_mut <- wt$r + dr

cat("THREE DIFFERENT HAMILTONIANS. Same structures, different springs.\n\n")
cat("H_wt   : springs l_wt        (the wild type)\n")
cat("H_mut  : springs l_wt + dl   (the LFENM mutant -- KEEPS the strain)\n")
cat("H_refit: springs l = d^e_mut (refit on the mutant structure)\n\n")

Vwt  <- venm(wt$r,  wt$pr, wt$k, wt$l)
Vmut_at_wt  <- venm(wt$r,  wt$pr, wt$k, l_mut)
Vmut_at_mut <- venm(r_mut, wt$pr, wt$k, l_mut)
rf <- build_relaxed(r_mut, wt$d_max, wt$k_val)
Vrefit_at_mut <- venm(rf$r, rf$pr, rf$k, rf$l)

cat(sprintf("V_wt   (r^e_wt)   = %10.6f   <- zero: wt is relaxed\n", Vwt))
cat(sprintf("V_mut  (r^e_wt)   = %10.6f   <- stress term, 1/2 sum k dl^2\n", Vmut_at_wt))
cat(sprintf("V_mut  (r^e_mut)  = %10.6f   <- stress MINUS relaxation. NOT zero.\n", Vmut_at_mut))
cat(sprintf("V_refit(r^e_mut)  = %10.6f   <- zero: refit springs are relaxed\n\n", Vrefit_at_mut))

cat("So the true energy change of the mutation is\n")
cat(sprintf("  dV = V_mut(r^e_mut) - V_wt(r^e_wt) = %.6f - %.6f = %.6f\n\n",
    Vmut_at_mut, Vwt, Vmut_at_mut - Vwt))
cat("and its two-term form (your sc_lfenm.md / wt_mut_energies.md):\n")
Vstress <- 0.5*sum(wt$k*dl^2)
Vrelax  <- 0.5*as.numeric(crossprod(dr, K %*% dr))
cat(sprintf("  1/2 sum k dl^2            = %10.6f\n", Vstress))
cat(sprintf("  1/2 (dr)' K (dr)          = %10.6f\n", Vrelax))
cat(sprintf("  difference                = %10.6f   (exact dV = %.6f, err %.1e)\n\n",
    Vstress - Vrelax, Vmut_at_mut - Vwt, (Vstress-Vrelax)-(Vmut_at_mut-Vwt)))
cat("The ERROR in the old text was to write V_mut(r^e_mut) = 0.\n")
cat("What is zero is V_refit(r^e_mut): a DIFFERENT Hamiltonian.\n")
