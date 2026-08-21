library(here)
## Is dV additive over sites? Catalogue says dV(s) = sum_j dV(0 -> j,s_j).
## Test against applying the same mutations TOGETHER and measuring exactly.
source(here("R", "catalog.R")); source(here("R", "fig_setup.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5); N <- length(r0)/3
K <- state_K(wt); C <- cmat_of(nma(K))

## apply a SET of site mutations at once, exactly
apply_set <- function(wt, sites, dls) {
  dl <- rep(0, nrow(wt$pr))
  for (q in seq_along(sites)) dl <- dl + dls[[q]]
  l2 <- wt$l + dl
  f <- force(wt$r, wt$pr, wt$k, l2); dr <- as.vector(C %*% f)
  st <- list(r=wt$r+dr, pr=wt$pr, k=wt$k, l=l2, frustrated=FALSE)
  rq <- state_req(st); st$r <- rq$r
  list(dV = venm(st$r,st$pr,st$k,l2) - venm(wt$r,wt$pr,wt$k,wt$l), dr = st$r - wt$r)
}
one <- function(wt, j, seed) {
  set.seed(seed); mask <- wt$pr[,1]==j | wt$pr[,2]==j
  dl <- rep(0, nrow(wt$pr)); dl[mask] <- rnorm(sum(mask),0,.3); dl
}
cat("Additivity of dV over simultaneous site mutations (2acy):\n\n")
cat(sprintf("%6s %12s %12s %10s %10s %8s\n","nsites","sum of parts","exact","abs err","rel err","shared"))
set.seed(1)
for (ns in c(2,3,5,10,20,40)) {
  errs <- rels <- shr <- numeric(6)
  for (rep in 1:6) {
    sites <- sample(N, ns)
    dls <- lapply(seq_along(sites), function(q) one(wt, sites[q], 1000*rep+q))
    ex <- apply_set(wt, sites, dls)
    parts <- sum(sapply(seq_along(sites), function(q) {
      l2 <- wt$l + dls[[q]]
      f <- force(wt$r,wt$pr,wt$k,l2); dr <- as.vector(C%*%f)
      st <- list(r=wt$r+dr,pr=wt$pr,k=wt$k,l=l2,frustrated=FALSE)
      st$r <- state_req(st)$r
      venm(st$r,st$pr,st$k,l2) - venm(wt$r,wt$pr,wt$k,wt$l) }))
    errs[rep] <- parts - ex$dV; rels[rep] <- (parts-ex$dV)/ex$dV
    ## how many contacts are shared between the mutated sites?
    m <- rowSums(sapply(dls, function(d) d != 0))
    shr[rep] <- sum(m > 1)
    if (rep==1 && ns %in% c(2,40))
      cat(sprintf("%6d %12.5f %12.5f %10.2e %9.2f%% %8d\n", ns, parts, ex$dV,
          parts-ex$dV, 100*(parts-ex$dV)/ex$dV, sum(m>1)))
  }
  if (!(ns %in% c(2,40)))
    cat(sprintf("%6d %12s %12s %10.2e %9.2f%% %8.1f\n", ns, "-","-",
        mean(errs), 100*mean(rels), mean(shr)))
}
