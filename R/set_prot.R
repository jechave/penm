#' Set up ENM from bio3d pdb object
#'
#' @param pdb a pdb object obtained using bio3d::read.pdb
#' @param node network nodes: "sc" or "ca"
#' @param model enm model: "anm", "ming_wall", "hnm", "hnm0", "pfanm", "reach"
#' @param d_max cutoff to define enm contacts
#' @param frustrated logical indicating whether to include frustrations in calculation of kmat
#'
#' @return a `prot` object (list containing xyz, site, etc.)
#' @export
#'
#' @examples
#'
set_prot <- function(pdb, node = "sc", model, d_max, frustrated) {

  prot <- nodes(pdb, node) # get xyz of nodes, pdb_site, and bfactor
  prot <- add_site_indexes(prot) # add site
  prot$enm_param <- lst(node, model, d_max, frustrated) # add enm parameters
  prot <- add_enm(prot, model, d_max, frustrated) # calculate enm and nma
  prot$energy <- enm_energy(prot) # calculate enm energy

  prot
}





nodes <- function(pdb, node) {
  if (node == "calpha" | node == "ca") {
    prot <- prot_ca(pdb)
  } else if (node == "side_chain" | node == "sc") {
    prot <- prot_sc(pdb)
  } else {
    stop("Error: node must be ca, calpha, sc, or side_chain")
  }
}

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


add_site_indexes <- function(prot) {
  prot <- prot
  prot$nsites <- length(prot$pdb_site) # add nsites
  prot$site <- seq(prot$nsites) # add site
  prot
}





