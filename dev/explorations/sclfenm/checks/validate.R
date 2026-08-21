library(here)
source(here("R", "sclfenm.R")); set.seed(1)
make_prot <- function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r0<-prot_2acy()$r
ok <- function(lbl,cond,extra="") cat(sprintf("[%s] %-52s %s\n", if(cond)"PASS" else "FAIL", lbl, extra))

cat("=== V1. energy is conserved across the rebuild ===\n")
set.seed(5); wt<-new_state(r0)
m_reb <- { set.seed(11); mutate(wt, 11, model="sclfenm") }
m_nor <- { set.seed(11); mutate(wt, 11, model="lfenm") }
e_reb <- state_vmin(m_reb$state); e_nor <- state_vmin(m_nor$state)
ok("same dV under lfenm and sclfenm", isTRUE(all.equal(m_reb$dV,m_nor$dV)),
   sprintf("dV=%.6f", m_reb$dV))
ok("total energy identical after rebuild", isTRUE(all.equal(e_reb,e_nor,tolerance=1e-10)),
   sprintf("%.6f vs %.6f", e_reb, e_nor))
ok("rebuilt network really is relaxed",
   isTRUE(all.equal(venm(m_reb$state$r,m_reb$state$pr,m_reb$state$k,m_reb$state$l),0)),
   "strain in network = 0, all energy in v_off")
cat("   v_off(reb) =", m_reb$state$v_off, "  v_off(nor) =", m_nor$state$v_off, "\n\n")

cat("=== V2. decomposition matches exact dV to O(dl^3) ===\n")
d <- m_reb$dV - (m_reb$V_stress - m_reb$V_relax)
ok("dV = V_stress - V_relax", abs(d) < 0.01*abs(m_reb$dV),
   sprintf("err=%.2e (%.3f%%)", d, 100*abs(d/m_reb$dV)))
cat("   V_stress =", m_reb$V_stress, " V_relax =", m_reb$V_relax, " dV =", m_reb$dV, "\n\n")

cat("=== V3. structure is the SAME under lfenm and sclfenm ===\n")
ok("mutant structures identical", isTRUE(all.equal(m_reb$state$r, m_nor$state$r)),
   sprintf("rmsd=%.2e", sqrt(mean((m_reb$state$r-m_nor$state$r)^2))))
cat("   (both relax with C_wt; the rebuild happens after)\n\n")

cat("=== V4. motion changes only under sclfenm ===\n")
mo_reb <- delta_motion(wt, m_reb$state); mo_nor <- delta_motion(wt, m_nor$state)
ok("lfenm: dTS is exactly zero", abs(mo_nor$dTS) < 1e-9, sprintf("dTS=%.3e", mo_nor$dTS))
ok("sclfenm: dTS is non-zero",  abs(mo_reb$dTS) > 1e-6, sprintf("dTS=%.6f", mo_reb$dTS))
cat("   rmsip lfenm =", round(mo_nor$rmsip,6), " sclfenm =", round(mo_reb$rmsip,6), "\n\n")

cat("=== V5. the founder has zero energy by construction ===\n")
ok("state_vmin(founder) == 0", isTRUE(all.equal(state_vmin(wt),0)))
