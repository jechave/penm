model_mutation <- function(xyz,pdb_site,sd_min,k_ij,site_to_mutate,d_max=12.5,dlij_max) {
  # Calculates the 'force' that models a mutation at a site
  nsites = length(xyz)/3
  site = seq(nsites) # use consecutive numbers for sites, regardless of pdb_site
  dim(xyz) <- c(3,nsites)
  # Sequence distance
  sdij = abs(pdb_site-pdb_site[site_to_mutate])
  # Distance
  dijV <- xyz - xyz[,site_to_mutate]
  dij <- sqrt(colSums(dijV^2))
  eij_versor <- t(t(dijV)/dij)
  f <- matrix(0,nrow=3,ncol=nsites)
  dv0 <- 0
  for (j in site) {
    if ( j != site_to_mutate & dij[j]<= d_max & sdij[j]>=sd_min) {
      # only perturb tertiary contacts
      dlij = runif(1,-dlij_max,dlij_max)
      fij <- k_ij[site_to_mutate,j]*dlij
      f[,j] = fij*eij_versor[,j]
      dv0 <- dv0 + .5* k_ij[site_to_mutate,j]*dlij^2
    }
  }
  f[,site_to_mutate] = -rowSums(f)
  dim(f) = c(length(f))
  ret = list("f"=f,"dv0" = dv0)
  ret
}
