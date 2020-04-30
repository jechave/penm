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

add_site_indexes <- function(prot) {
  prot <- prot
  prot$nsites <- length(prot$pdb_site) # add nsites
  prot$site <- seq(prot$nsites) # add site
  prot
}



