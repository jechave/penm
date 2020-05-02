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
  prot$enm_param <- lst(node, model, d_max, frustrated) # add enm parameters
  prot <- add_enm(prot, model, d_max, frustrated) # calculate enm and nma

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
  nsites <- length( r$site )
  site <- seq( length( r$site ) )

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


  result <- lst( nsites,
                 site,
                 pdb_site,
                 bfactor,
                 xyz )
  result
}

prot_ca <- function(pdb) {
  sel <- atom.select(pdb, elety = "CA") # select CA

  nsites <- length(sel$atom)
  site <- seq(length(sel$atom))
  pdb_site <- pdb$atom$resno[sel$atom]
  bfactor <- pdb$atom$b[sel$atom]
  xyz <- pdb$xyz[sel$xyz]

  result <- lst( nsites,
                 site,
                 pdb_site,
                 bfactor,
                 xyz )
  result
}


site <- function(prot) {
  seq( length(prot$pdb_site) )
}

nsites <- function(prot) {
  length(prot$pdb_site)
}





