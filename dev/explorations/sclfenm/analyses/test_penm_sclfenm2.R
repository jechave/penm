## find the package root by walking up, so this survives being moved
.penm_root <- local({ d <- normalizePath(".")
  while (!file.exists(file.path(d, "DESCRIPTION")) && dirname(d) != d) d <- dirname(d); d })
suppressPackageStartupMessages(devtools::load_all(.penm_root, quiet = TRUE))
wt <- set_enm(pdb_2acy_A, node="ca", model="anm", d_max=10.5, frustrated=FALSE)
enm_f <- function(p) { g<-get_graph(p); R<-matrix(as.vector(get_xyz(p)),nrow=3)
  D<-R[,g$j,drop=FALSE]-R[,g$i,drop=FALSE]; d<-sqrt(colSums(D^2)); E<-sweep(D,2,d,"/")
  a<-g$kij*(d-g$lij); f<-matrix(0,3,ncol(R))
  for(q in seq_len(nrow(g))){f[,g$i[q]]<-f[,g$i[q]]+a[q]*E[,q]; f[,g$j[q]]<-f[,g$j[q]]-a[q]*E[,q]}
  sqrt(sum(f^2)) }

cat("Why is the sclfenm mutant so far from its own minimum?\n")
cat("Hypothesis: the structure is computed by ONE linear step with C_wt, then\n")
cat("the graph is rebuilt around it -- but the structure is never re-relaxed\n")
cat("against the NEW graph. So r is a minimum of neither graph.\n\n")
cat(sprintf("%6s %14s %14s %12s\n","site","|F| lfenm","|F| sclfenm","ratio"))
set.seed(2)
for (j in sample(get_nsites(wt), 6)) {
  a <- enm_f(get_mutant_site(wt,j,1,"lfenm",.3,1L,1L))
  b <- enm_f(get_mutant_site(wt,j,1,"sclfenm",.3,1L,1L))
  cat(sprintf("%6d %14.4e %14.4e %12.1f\n", j, a, b, b/a))
}
cat("\nWhere does the extra force come from? Compare the mutant's graph lij\n")
cat("against the distances in its own structure:\n")
mu <- get_mutant_site(wt,40,1,"sclfenm",.3,1L,1L); g <- get_graph(mu)
new_edges <- !(g$edge %in% get_graph(wt)$edge)
cat(sprintf("  edges kept from wt : %d, max|d-l| = %.4f\n",
    sum(!new_edges), max(abs(g$dij-g$lij)[!new_edges])))
cat(sprintf("  edges NEW in mutant: %d, max|d-l| = %.4f  <- relaxed by construction\n",
    sum(new_edges), if(any(new_edges)) max(abs(g$dij-g$lij)[new_edges]) else 0))
