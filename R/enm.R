
# Create and set prot object ----------------------------------------------


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
set_enm <- function(pdb, node, model, d_max, frustrated) {
  prot <- create_enm() %>%
    set_enm_param(node = node, model = model, d_max = d_max, frustrated = frustrated) %>%
    set_enm_nodes(pdb = pdb) %>%
    set_enm_graph() %>%
    set_enm_eij() %>%
    set_enm_kmat() %>%
    set_enm_nma()

  prot
}



# Set prot components ------------------------------------------------

#' Create an empty prot object
#'
create_enm <- function() {
  prot <- lst(param = NA, nodes = NA, graph = NA, eij = NA, kmat = NA, nma = NA)
  class(prot) <- "prot"
  prot
}

#' Set param of prot object
#'
set_enm_param <- function(prot, node, model, d_max, frustrated) {
  prot$param <- lst(node, model, d_max, frustrated)
  prot
}


#' Set nodes of prot object
#'
set_enm_nodes <- function(prot, pdb) {
  prot$nodes <- calculate_enm_nodes(pdb, get_enm_node(prot))
  return(prot)
}


#' Set graph of prot object
#'
set_enm_graph <- function(prot) {
  prot$graph <- calculate_enm_graph(get_xyz(prot), get_pdb_site(prot), get_enm_model(prot), get_d_max(prot))
  prot
}

#' Set eij unit vectors of prot object
#'
set_enm_eij <- function(prot) {
  prot$eij <- calculate_enm_eij(get_xyz(prot), get_graph(prot)$i, get_graph(prot)$j)
  prot
}

#' Set enm's kmat of prot object
#'
set_enm_kmat <- function(prot) {
  prot$kmat <- calculate_enm_kmat(get_graph(prot), get_eij(prot), get_nsites(prot), get_frustrated(prot))
  prot
}

#' Set normal-mode-analysis component of prot object
set_enm_nma <- function(prot) {
  prot$nma <- calculate_enm_nma(get_kmat(prot))
  prot
}


# Calculate prot components -----------------------------------------------


#' Calculate nodes of prot object
#'
calculate_enm_nodes <- function(pdb, node) {
  if (node == "calpha" | node == "ca") {
    nodes <- prot_ca(pdb)
    return(nodes)
  }
  if (node == "side_chain" | node == "sc") {
    nodes <- prot_sc(pdb)
    return(nodes)
  }
  stop("Error: node must be ca, calpha, sc, or side_chain")
}

#' set side-chain nodes
#'
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


  lst( nsites, site, pdb_site, bfactor, xyz )
}

#' set CA nodes
#'
prot_ca <- function(pdb) {
  sel <- atom.select(pdb, elety = "CA") # select CA

  nsites <- length(sel$atom)
  site <- seq(length(sel$atom))
  pdb_site <- pdb$atom$resno[sel$atom]
  bfactor <- pdb$atom$b[sel$atom]
  xyz <- pdb$xyz[sel$xyz]

  lst( nsites, site, pdb_site, bfactor, xyz )
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
#' if (FALSE)
#'  calculate_enm_graph(xyz, pdb_site, model, d_max)
#'
#'@family enm builders
calculate_enm_graph <- function(xyz, pdb_site, model, d_max,...) {
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
#'
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
#'
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
calculate_enm_eij <- function(xyz,i,j) {

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
#' calculate_enm_kmat(graph, eij)
#'
#' @export
#'
#'@family enm builders
#'
calculate_enm_kmat <- function(graph, eij, nsites, frustrated) {
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


#' Perform Normal Mode Analysis
#'
#' Given an enm `kmat`, perform NMA
#'
#' @param kmat The K matrix to diagonalize
#' @param TOL=1.e-5 A small value, eigenvectors with eigenvalues larger than `TOL` are discarded
#'
#' @return A list with elements \code{mode}, \code{evalue}, \code{cmat}, and \code{umat}
#'
#' @examples
#'
#' enm_anm(enm, 1.e-10)
#'
#' @export
#'
#'@family enm builders
#'
calculate_enm_nma <- function(kmat, TOL = 1.e-5) {
  # Given an enm object for which kmat has already been defined, perform NMA
  # It returns a list containing mode, evalue, cmat, umat
  # Diagonalize
  eig <- eigen(kmat, symmetric = TRUE)
  evalue <- eig$values
  umat <- eig$vectors
  modes <- evalue > TOL
  evalue <- evalue[modes]
  umat  <- umat[, modes]

  nmodes <- sum(modes)
  mode <- seq.int(from = nmodes, to = 1, by = -1)
  evalue <- evalue[mode]
  umat <- umat[, mode]
  mode <- mode[mode]


  cmat <-  umat %*% ((1 / evalue) * t(umat))

  nma <- list(
    mode = mode,
    evalue = evalue,
    cmat = cmat,
    umat = umat
  )
  nma
}



