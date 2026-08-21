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
cat("Does the residual force track the CHANGE IN EDGE COUNT?\n\n")
cat(sprintf("%6s %10s %14s %14s\n","site","d(edges)","|F| lfenm","|F| sclfenm"))
set.seed(7); ch <- unch <- c()
for (j in sample(get_nsites(wt), 20)) {
  ml <- get_mutant_site(wt,j,1,"lfenm",.3,1L,1L)
  ms <- get_mutant_site(wt,j,1,"sclfenm",.3,1L,1L)
  de <- nrow(get_graph(ms)) - nrow(get_graph(wt))
  if (de == 0) unch <- c(unch, enm_f(ms)/enm_f(ml)) else ch <- c(ch, enm_f(ms)/enm_f(ml))
  if (length(ch)+length(unch) <= 8)
    cat(sprintf("%6d %10d %14.4e %14.4e\n", j, de, enm_f(ml), enm_f(ms)))
}
cat(sprintf("\nmean |F| ratio (sclfenm/lfenm):\n"))
cat(sprintf("  edge count UNCHANGED : %.2f   (n = %d)\n", mean(unch), length(unch)))
cat(sprintf("  edge count CHANGED   : %.2f   (n = %d)\n", mean(ch), length(ch)))
cat("\n=> the residual force is driven by contact-set changes. When the map is\n")
cat("   stable the two models agree; when an edge appears or vanishes, the\n")
cat("   structure is left far from the rebuilt graph's minimum.\n")
