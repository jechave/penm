

prot_sc <- function(pdb) {
  # calculate xyz coordinates of ca, cb, com, qm (qb by Micheletti), ql (qb by Levitt)
  r = residue.coordinates(pdb)
  b = residue.bfactors(pdb)

  pdb_site <- r$site

  xyz <-  r$com.xyz
  xyz_na <- is.na(xyz)
  xyz[xyz_na] <- r$m.xyz[xyz_na] #when center of mass is NA, use Michelettis approximation to Cbeta coordinates
  xyz_na <- is.na(xyz)
  xyz[xyz_na] <- r$ca.xyz[xyz_na] # use CA coordinates if neither com or micheletti are defined

  b.c <- b$b.c
  b.m <- b$b.m
  b.a <- b$b.a

  bfactor <- b.c # center-of-mass bfactors
  bfactor[is.na(bfactor)] <- b.m[is.na(bfactor)] # Use Micheletti's bfactor when com bfactor undefined'
  bfactor[is.na(bfactor)] <- b.a[is.na(bfactor)] # Use CA bfactor otherwise


  lst(xyz, pdb_site, bfactor)
}

prot_ca <- function(pdb) {
  # select CA
  sel <- atom.select(pdb, elety = "CA")
  # prepare output
  bfactor <- as.numeric(pdb$atom[sel$atom, c("b")])
  pdb_site <- as.numeric(pdb$atom[sel$atom, c("resno")])
  nsites <- length(pdb_site)
  xyz.calpha <-
    matrix(pdb$xyz[sel$xyz],
           ncol = nsites,
           nrow = 3,
           byrow = F)

  output <-
    list(
      "xyz" = as.vector(xyz.calpha),
      "pdb_site" = pdb_site,
      "bfactor" = bfactor
    )
  output
}
