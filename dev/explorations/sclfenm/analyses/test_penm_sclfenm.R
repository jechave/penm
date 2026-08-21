## What does this exploration say about penm's own sclfenm/frustration handling?
## find the package root by walking up, so this survives being moved
.penm_root <- local({ d <- normalizePath(".")
  while (!file.exists(file.path(d, "DESCRIPTION")) && dirname(d) != d) d <- dirname(d); d })
suppressPackageStartupMessages(devtools::load_all(.penm_root, quiet = TRUE))
wt <- set_enm(pdb_2acy_A, node="ca", model="anm", d_max=10.5, frustrated=FALSE)
N <- get_nsites(wt)

cat("=== ISSUE 1: mutate_graph()'s lij copy ===\n")
cat("The line:  g2[g2$edge %in% g1$edge, 'lij'] <- g1[g1$edge %in% g2$edge, 'lij']\n")
cat("carries WARNING: 'true only if edges are ordered'. Test whether the two\n")
cat("selections actually align when the contact set changes.\n\n")
mu <- get_mutant_site(wt, 40, 1, "sclfenm", .3, 1L, 1L)
g1 <- get_graph(wt); g2 <- get_graph(mu)
cat(sprintf("  wt edges %d ; mutant edges %d ; changed? %s\n",
    nrow(g1), nrow(g2), nrow(g1)!=nrow(g2)))
lhs <- which(g2$edge %in% g1$edge); rhs <- which(g1$edge %in% g2$edge)
cat(sprintf("  LHS selects %d rows, RHS selects %d rows\n", length(lhs), length(rhs)))
cat(sprintf("  do the selected edge LABELS match pairwise? %s\n",
    identical(g2$edge[lhs], g1$edge[rhs])))
cat("  -> if TRUE the copy is correct here; the warning is about the general case\n\n")

cat("=== ISSUE 2: is the mutant's own dij consistent with its xyz? ===\n")
d_true <- dij_edge(matrix(as.vector(get_xyz(mu)),nrow=3), g2$i, g2$j)
cat(sprintf("  max|graph$dij - dij(xyz)| = %.3e\n", max(abs(g2$dij - d_true))))

cat("\n=== ISSUE 3: is eij refreshed? (it is used to build kmat) ===\n")
e_stored <- get_eij(mu)
e_true <- calculate_enm_eij(matrix(as.vector(get_xyz(mu)),nrow=3), g2$i, g2$j)
cat(sprintf("  max|eij stored - eij(xyz)| = %.3e\n", max(abs(e_stored - e_true))))
cat("  (sclfenm calls set_enm_eij(), so this SHOULD be 0)\n")

cat("\n=== ISSUE 4: does sclfenm leave the mutant at its own minimum? ===\n")
enm_f <- function(p) { g<-get_graph(p); R<-matrix(as.vector(get_xyz(p)),nrow=3)
  D<-R[,g$j,drop=FALSE]-R[,g$i,drop=FALSE]; d<-sqrt(colSums(D^2)); E<-sweep(D,2,d,"/")
  a<-g$kij*(d-g$lij); f<-matrix(0,3,ncol(R))
  for(q in seq_len(nrow(g))){f[,g$i[q]]<-f[,g$i[q]]+a[q]*E[,q]; f[,g$j[q]]<-f[,g$j[q]]-a[q]*E[,q]}
  sqrt(sum(f^2)) }
cat(sprintf("  |F| of the sclfenm mutant on its OWN graph = %.4e\n", enm_f(mu)))
mu_l <- get_mutant_site(wt, 40, 1, "lfenm", .3, 1L, 1L)
cat(sprintf("  |F| of the lfenm   mutant on its own graph = %.4e\n", enm_f(mu_l)))
