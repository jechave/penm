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
  prot$enm <- enm_from_xyz(prot$xyz, prot$pdb_site, model, d_max, frustrated)

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









#' Set up ENM model
#'
#' Set up ENM model from node coordinates `xyz` and residue pdb numbers `pdb_site`.
#'
#' @param xyz The \code{3 x N} matrix containing the Cartesian coordinates of ENM nodes
#' @param pdb_site The pdb numbering of network nodes
#' @param model The ENM model variant (gnm, anm, hnm0, hnm, ming_wall, reach, pfgnm)
#' @param d_max  The cut-off used to define contacts in some models
#' @param ... Any other parameter
#'
#' @return A list containing the ENM graph `graph`, the matrix of versors `eij`, and the ENM K matrix `kmat`
#'
#' @export
#'
#' @examples
#'
#' @family enm builders
enm_from_xyz <- function(xyz, pdb_site, model,  d_max, frustrated,...) {

  # Calculate (relaxed) enm graph from xyz
  # Returns list (graph, eij, kmat)

  graph <- enm_graph_xyz(xyz, pdb_site, model,  d_max)

  # calculate eij
  eij <- eij_edge(xyz, graph$i, graph$j)

  # calculate kmat
  kmat <- kmat_graph(graph, eij, nsites = length(pdb_site), frustrated)

  # diagonalise kmat
  nma <- enm_nma(kmat)

  # add (mode, evalue, umat, cmat)

  # return enm object
  result <- lst(graph = graph,
                eij = eij,
                kmat = kmat,
                mode = nma$mode,
                evalue = nma$evalue,
                cmat = nma$cmat,
                umat = nma$umat)

  result
}

#' Calculate ENM graph
#'
#' Calculates graph representation of Elastic Network Model (ENM), the typical relaxed case (lij = dij)
#'
#' @param xyz matrix of size \code{c(3,N)} containing each column the \code{x, y, z} coordinates of each of N nodes
#' @param pdb_site integer vector of size N containing the number of each node (pdb residue number)
#' @param model  character variable specifying the ENM model variant, default is \code{"gnm"}, options:
#'     \code{gnm, anm, ming_wall, hnm, hnm0, pfgnm, reach}.
#' @param d_max distance-cutoff to define network contacts
#' @return a tibble that contains the graph representation of the network
#'
#' @export
#'
#' @examples
#' enm_graph_xyz(xyz, pdb_site)
#'
#'@family enm builders
enm_graph_xyz <- function(xyz,
           pdb_site,
           model,
           d_max,
           ...) {
    # Calculate (relaxed) enm graph from xyz
    # Returns graph for the relaxed case

    # put xyz in the right format and check size
    xyz <- my_as_xyz(xyz)
    nsites <- length(pdb_site)
    stopifnot(ncol(xyz) == nsites)

    # set function to calculate i-j spring constants
    kij_fun <- match.fun(paste0("kij_", model))
    kij_par <- lst( d_max = d_max)

    site <- seq(nsites)
    # calculate graph
    graph <- as_tibble(expand.grid(i = site, j = site)) %>%
      filter(j > i) %>%
      arrange(i, j) %>%
      mutate(dij = dij_edge(xyz, i, j)) %>%
      mutate(sdij = sdij_edge(pdb_site, i, j)) %>%
      filter(dij <= d_max | sdij == 1) %>%
      mutate(lij = dij)

    graph$kij <- do.call(kij_fun,
                         c(lst(
                           dij = graph$dij, sdij = graph$sdij
                         ), kij_par))

    graph <- graph %>%
      mutate(edge = paste(i, j, sep = "-"),
             lij = dij) %>%
      mutate(v0ij = 0) %>%
      dplyr::select(edge, i, j, v0ij, sdij, lij, kij, dij)

    graph
  }

#' Calculate distance of edges
dij_edge <- function(xyz,i,j) {
  stopifnot(length(i) == length(j))
  xyz <- my_as_xyz(xyz)
  dij <- rep(NA,length(i))
  for (k in seq(length(i)))  {
    ik <- i[k]
    jk <- j[k]
    rij <- xyz[,jk] - xyz[,ik]
    dij[k] <- sqrt(sum(rij^2))
  }
  dij
}

#' Calculate edge sequence distance
sdij_edge <- function(pdb_site,i,j) {
  # sequence distance
  stopifnot(length(i) == length(j))
  n_edges <- length(i)
  sdij <- abs(pdb_site[j] - pdb_site[i])
  sdij

}


#' Calculate unit vectors of edges
#'
#' @param i,j integer vectors of nodes connected in each edge
#' @param xyz vector of xyz coordinates
#' @return matrix with n_edge rows and 3 columns (x, y, z)
#'
eij_edge <- function(xyz,i,j) {

  stopifnot(length(i) == length(j))
  n_edges <- length(i)
  xyz <- my_as_xyz(xyz)
  # eij <- tibble(eij_x = rep(NA,n_edges),
  #               eij_y = rep(NA,n_edges),
  #               eij_z = rep(NA,n_edges))

  eij <- matrix(NA,n_edges,3)
  for (k in seq(n_edges)) {
    ik <- i[k]
    jk <- j[k]
    rij <- xyz[,jk] - xyz[,ik]
    eij[k,] <- rij/sqrt(sum(rij^2))
  }
  eij
}




#' Calculate kmat given the ENM graph
#'
#' @param graph A tibble representing the ENM graph (with edge information, especially \code{kij}
#' @param eij A matrix of size \code{n_edges x 3} of \code{eij} versors directed along ENM contacts
#' @param nsites The number of nodes of the ENM network
#' @param frustrated Logical indicating  wether to add frustration or not before calculating \code{kmat}
#'
#' @return The \code{3 nsites x 3 nsites} stiffness matrix of the ENM
#'
#' @examples
#' kmat_graph(graph, eij)
#'
#' @export
#'
#'@family enm builders
#'
kmat_graph <- function(graph, eij, nsites, frustrated) {
  stopifnot(max(graph$i, graph$j) <= nsites,
            nrow(graph) == nrow(eij))
  kmat <- array(0, dim = c(3, nsites, 3, nsites))
  for (edge in seq(nrow(graph))) {
    i <- graph$i[[edge]]
    j <- graph$j[[edge]]
    kij <- graph$kij[[edge]]
    if (frustrated) {
      gij <- graph$lij[[edge]]/graph$dij[[edge]] - 1
    } else {
      gij <- 0
    }
    eij_v <- eij[edge,]
    eij_mat <- tcrossprod(eij_v, eij_v)
    kij_mat <- -kij * (eij_mat + gij * (eij_mat - diag(3)))
    kmat[, j, , i] <- kmat[, i, , j] <- kij_mat
  }
  for (i in seq(nsites)) {
    kmat[, i, , i] <- -apply(kmat[, i, , -i], c(1, 2), sum)
  }

  dim(kmat) <- c(3 * nsites, 3 * nsites)
  kmat
}


