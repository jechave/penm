library(here)
source(here("R", "sclfenm.R")); source(here("R", "fig_setup.R"))
## proper minimiser: iterate LRA until the force vanishes
minimise <- function(st, ...) state_req(st)   # use the model's own solver
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5); N <- length(r0)/3
cat("How wrong is V(st$r) compared with the true minimum V(r_eq)?\n\n")
cat(sprintf("%-10s %6s %12s %12s %12s %10s\n",
    "model","nmut","V(st$r)","V(r_eq)","difference","rel.err"))
for (md in c("keepnet","lfenm")) {
  st <- wt; set.seed(3)
  for (n in 1:10) {
    st <- mutate(st, sample(N,1), sigma=.3, model=md)$state
    if (n %in% c(1,3,5,10)) {
      m <- minimise(st)
      Vr  <- venm(st$r, st$pr, st$k, st$l) + st$v_off
      Veq <- venm(m$r,  st$pr, st$k, st$l) + st$v_off
      cat(sprintf("%-10s %6d %12.5f %12.5f %12.3e %9.2f%%\n",
          md, n, Vr, Veq, Vr-Veq, 100*(Vr-Veq)/Veq))
    }
  }
}
cat("\nsclfenm (r is exact by construction):\n")
st <- wt; set.seed(3)
for (n in 1:10) {
  st <- mutate(st, sample(N,1), sigma=.3, model="sclfenm")$state
  if (n==10) { m<-minimise(st)
    cat(sprintf("%-10s %6d %12.5f %12.5f %12.3e\n","sclfenm",n,
      venm(st$r,st$pr,st$k,st$l)+st$v_off, venm(m$r,st$pr,st$k,st$l)+st$v_off,
      venm(st$r,st$pr,st$k,st$l)-venm(m$r,st$pr,st$k,st$l))) }
}
