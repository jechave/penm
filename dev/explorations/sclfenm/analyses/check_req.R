library(here)
source(here("R", "sclfenm.R")); source(here("R", "fig_setup.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5); N <- length(r0)/3
cat("Is st$r the equilibrium conformation of st's own network?\n")
cat("Test: residual force |F(st$r)| under st's own (pr,k,l).\n\n")
show <- function(tag, st) {
  f <- force(st$r, st$pr, st$k, st$l)
  strain <- 0.5*sum(st$k*(dists(st$r,st$pr)-st$l)^2)
  cat(sprintf("%-28s |F| = %10.3e   strain at st$r = %9.5f   v_off = %8.5f\n",
              tag, sqrt(sum(f^2)), strain, st$v_off))
}
show("founder", wt)
set.seed(11)
for (md in c("lfenm","keepnet","sclfenm")) {
  mu <- mutate(wt, 11, sigma=.3, model=md)
  show(paste0("after 1 mutation: ", md), mu$state)
}
cat("\nSame, after 5 mutations:\n")
for (md in c("keepnet","sclfenm")) {
  st <- wt; set.seed(3)
  for (i in 1:5) st <- mutate(st, sample(N,1), sigma=.3, model=md)$state
  show(paste0("5 mutations: ", md), st)
}
