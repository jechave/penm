library(here)
## Background-independence, both models, same protocol, several trajectory lengths.
source(here("R", "frustration.R"))
wt <- founder(); N <- get_nsites(wt)
effect <- function(prot, dl) {                 # apply a GIVEN dl to a background
  K <- hess(prot, FALSE); s <- spec(K); C <- s$vec %*% ((1/s$val)*t(s$vec))
  f <- calculate_force(prot, dl); dr <- as.vector(C %*% f)
  g <- get_graph(prot); l2 <- g$lij + dl
  xyz <- matrix(as.vector(get_xyz(prot)),nrow=3) + matrix(dr,nrow=3)
  d2 <- dij_edge(xyz, g$i, g$j)
  list(dV = 0.5*sum(g$kij*(d2-l2)^2) - 0.5*sum(g$kij*(g$dij-g$lij)^2),
       dr = dr)
}
j <- 40; set.seed(99); mask <- get_graph(wt)$i==j | get_graph(wt)$j==j
dls <- lapply(1:8, function(m){d<-rep(0,nrow(get_graph(wt))); d[mask]<-rnorm(sum(mask),0,.3); d})
b0 <- lapply(dls, function(d) effect(wt, d))
aV <- sapply(b0,`[[`,"dV"); aR <- sapply(b0, function(x) sqrt(sum(x$dr^2)))
cat("SC-LFENM cannot be catalogued at all: the rebuild changes the edge count,\n")
cat("so a dl built on the founder no longer matches the mutant graph.\n")
mu <- wt; set.seed(7)
for (n in 1:5) mu <- get_mutant_site(mu, sample(setdiff(1:N,j),1), n, "sclfenm", .3, 1L, 1L)
cat(sprintf("  founder edges = %d ; after 5 sclfenm steps = %d\n\n",
    nrow(get_graph(wt)), nrow(get_graph(mu))))
cat(sprintf("%-9s %6s %10s %12s %12s\n","model","nmut","cor(dV)","rel err dV","rel err |dr|"))
for (md in c("lfenm")) {
  mu <- wt; set.seed(7)
  for (n in 1:40) {
    mu <- get_mutant_site(mu, sample(setdiff(1:N,j),1), n, md, .3, 1L, 1L)
    if (n %in% c(5,10,20,40)) {
      b1 <- lapply(dls, function(d) effect(mu, d))
      bV <- sapply(b1,`[[`,"dV"); bR <- sapply(b1, function(x) sqrt(sum(x$dr^2)))
      cat(sprintf("%-9s %6d %10.4f %11.1f%% %11.2f%%\n", md, n, cor(aV,bV),
          100*mean(abs(bV-aV)/abs(aV)), 100*mean(abs(bR-aR)/aR)))
    }
  }
}
