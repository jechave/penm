library(here)
## Do my reimplementation (enm_core.R) and penm agree on the same problem?
## Sections 1-7 use mine; 8-10 use penm's. They should give the same numbers.
## find the package root by walking up, so this survives being moved
.penm_root <- local({ d <- normalizePath(".")
  while (!file.exists(file.path(d, "DESCRIPTION")) && dirname(d) != d) d <- dirname(d); d })
suppressPackageStartupMessages(devtools::load_all(.penm_root, quiet = TRUE))
source(here("R", "enm_core.R"))

## penm's wild type (ANM, d_max 10.5) vs mine built on the SAME coordinates
p  <- set_enm(pdb_2acy_A, node="ca", model="anm", d_max=10.5, frustrated=FALSE)
xyz <- as.vector(get_xyz(p))
mine <- build_relaxed(xyz, d_max=10.5, k_val=1)

cat("network:\n")
cat(sprintf("  edges   penm %4d   mine %4d   same? %s\n",
    nrow(get_graph(p)), nrow(mine$pr),
    nrow(get_graph(p))==nrow(mine$pr)))
gp <- get_graph(p)
same_pairs <- all(sort(paste(gp$i,gp$j)) == sort(paste(mine$pr[,1],mine$pr[,2])))
cat(sprintf("  identical contact set?  %s\n", same_pairs))
cat(sprintf("  max|l_penm - l_mine|    %.3e\n",
    max(abs(sort(gp$lij) - sort(mine$l)))))

## Hessian
Kp <- get_kmat(p)
Km <- hessian(mine$r, mine$pr, mine$k, mine$l, frustrated=FALSE)
cat(sprintf("\nHessian: max|K_penm - K_mine| = %.3e   (scale %.2f)\n",
    max(abs(Kp - Km)), max(abs(Kp))))

## spectrum
ep <- get_evalue(p); em <- nma(Km)$evalue
n <- min(length(ep), length(em))
cat(sprintf("modes: penm %d, mine %d ; max|eval diff| = %.3e\n",
    length(ep), length(em), max(abs(sort(ep)[1:n] - sort(em)[1:n]))))

## a mutation, same dl through both paths
set.seed(1)
mask_p <- gp$i==40 | gp$j==40
dl <- rep(0, nrow(gp)); dl[mask_p] <- rnorm(sum(mask_p), 0, .3)
f_p <- calculate_force(p, dl)
f_m <- force(mine$r, mine$pr, mine$k, mine$l + dl)
cat(sprintf("\nforce from the same dl: max|f_penm - f_mine| = %.3e (|f| = %.3f)\n",
    max(abs(f_p - f_m)), sqrt(sum(f_p^2))))
