library(here)
## The catalogue's ONE assumption: a site's states have background-independent
## effects. Catalogue site j against the founder; then re-measure the SAME
## perturbations dl against an evolved background and compare.
## (The difference identity dV(m0->m1)=dV(0->m1)-dV(0->m0) is exact by
##  construction -- energy is a state function -- so it is not tested here.)
source(here("R", "catalog.R")); source(here("R", "fig_setup.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5); N <- length(r0)/3

## measure the effect of a GIVEN dl on a GIVEN background state
effect <- function(st, dl) {
  K <- state_K(st); C <- cmat_of(nma(K))
  l2 <- st$l + dl
  f <- force(st$r, st$pr, st$k, l2); dr <- as.vector(C %*% f)
  s2 <- list(r=st$r+dr, pr=st$pr, k=st$k, l=l2, frustrated=st$frustrated)
  s2$r <- state_req(s2)$r
  list(dV = venm(s2$r,s2$pr,s2$k,l2) - venm(st$r,st$pr,st$k,st$l),
       dr2 = sum((s2$r-st$r)^2))
}
## a fixed set of perturbations for site j (the "states" of that site)
states_of <- function(st, j, M=6, seed=1) {
  set.seed(seed); mask <- st$pr[,1]==j | st$pr[,2]==j
  lapply(seq_len(M), function(m){ d<-rep(0,nrow(st$pr))
    d[mask]<-rnorm(sum(mask),0,.3); d })
}
j <- 40
dls <- states_of(wt, j, M=6, seed=99)
base <- lapply(dls, function(d) effect(wt, d))

cat("Site", j, ": 6 states catalogued on the FOUNDER, then re-measured on\n")
cat("backgrounds evolved 5/10/20/40 substitutions away (sclfenm).\n\n")
cat(sprintf("%-10s %8s %10s %10s %10s\n","background","steps","cor(dV)","mean|ddV|","rel err"))
cat(sprintf("%-10s %8d %10.4f %10.5f %9.1f%%\n","founder",0,1,0,0))
st <- wt; set.seed(7)
for (n in 1:40) {
  st <- mutate(st, sample(setdiff(1:N,j),1), sigma=.3, model="sclfenm")$state
  if (n %in% c(5,10,20,40)) {
    now <- lapply(dls, function(d) effect(st, d))
    a <- sapply(base,`[[`,"dV"); b <- sapply(now,`[[`,"dV")
    cat(sprintf("%-10s %8d %10.4f %10.5f %9.1f%%\n","evolved",n,
        cor(a,b), mean(abs(b-a)), 100*mean(abs(b-a)/abs(a))))
  }
}
