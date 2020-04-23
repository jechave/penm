#' Add \code{enm} object to \code{prot} object
#'
#' @param prot A protein object that must contain \code{xyz} and \code{pdb_site} elements, and active_site indexes \code{ind_active} (which may be NA)
#' @param model the enm model
#' @param v0 the energy at the minimum (for the relaxed model)
#' @param d_max a cdistance cut-off needed by some models to define sites in contact
#'
#' @return A protein object equal to input with enm added, where enm is a list containing \code{graph, eij, kmat, mode, evalue, cmat, umat}
#'  and also \code{cmat_active, kmat_active}, which are NA if \code{ind_active} is NA.
#'
#' @export
#'
#' @family enm builders
#'
#' @examples
enm_add <- function(prot, model, v0,  d_max,...)  {
  stopifnot(!is.null(prot$ind_active)) # stop if ind_active undefined (but not if it's NA)
  stopifnot(is.null(prot$enm$umat)) # it adds nma only if not already defined
  prot$enm <- enm_set(prot, model, v0, d_max) #sets graph and eij
  nma <- enm_nma(prot$enm)
  prot$enm <- c(prot$enm, nma) #append (mode, evalue, umat, cmat)

  if (anyNA(prot$ind_active)) { # these matrices are only defined if ind_active is defined
    prot$enm$cmat_active <- NA
    prot$enm$kmat_active <- NA
  } else {
    prot$enm$cmat_active <- prot$enm$cmat[prot$ind_active, prot$ind_active]
    prot$enm$kmat_active <- solve(prot$enm$cmat_active)
  }

  prot
}

#' @rdname enm_add
#' @export
add_enm <- enm_add

#' Set initial \code{enm} object
#'
#' @param prot A protein object that must contain \code{xyz} and \code{pdb_site} elements, and active_site indexes \code{ind_active} (which may be NA)
#' @param model the enm model
#' @param v0 the energy at the minimum (for the relaxed model)
#' @param d_max a cdistance cut-off needed by some models to define sites in contact
#'
#' @return A list containing the ENM graph `graph`, the matrix of versors `eij`, and the ENM K matrix `kmat`
#' @export
#'
#' @family enm builders
#'
#' @examples
enm_set <- function(prot, model, v0,  d_max,...) {
  vars = lst(xyz = prot$xyz, pdb_site = prot$pdb_site)
  param <- lst(model = model,  v0 = v0, d_max = d_max)
  args = c(vars,param)
  enm <- do.call(enm_set_xyz, args)
  enm
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
enm_set_xyz <- function(xyz, pdb_site, model,  d_max,...) {

  # Calculate (relaxed) enm graph from xyz
  # Returns list (graph, eij, kmat)

  graph <- enm_graph_xyz(xyz, pdb_site, model,  d_max)

  # calculate eij
  eij <- eij_edge(xyz, graph$i, graph$j)

  # calculate kmat
  kmat <- kmat_graph(graph, eij, nsites = length(pdb_site))

  # return enm object
  lst(graph = graph, eij = eij, kmat = kmat)
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
enm_graph_xyz <-
  function(xyz,
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
      filter(dij <= d_max) %>%
      mutate(sdij = sdij_edge(pdb_site, i, j),
             lij = dij)
    graph$kij <- do.call(kij_fun,
                         c(lst(
                           dij = graph$dij, sdij = graph$sdij
                         ), kij_par))

    graph <- graph %>%
      mutate(edge = paste(i, j, sep = "-"),
             lij = dij) %>%
      dplyr::select(edge, i, j, sdij, lij, kij, dij)

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
#' @param add_frust=FALSE Wether to add frustration or not before calculating \code{kmat}
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
kmat_graph <- function(graph, eij, nsites, add_frust = FALSE) {
  stopifnot(max(graph$i, graph$j) <= nsites,
            nrow(graph) == nrow(eij))
  kmat <- array(0, dim = c(3, nsites, 3, nsites))
  for (edge in seq(nrow(graph))) {
    i <- graph$i[[edge]]
    j <- graph$j[[edge]]
    kij <- graph$kij[[edge]]
    if (add_frust) {
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


#' Add v0 to ENM
#'
#' Add \code{v0ij = v0 / nedges} to each edge of graph of enm object
#' @param prot A protein object
#' @param v0 A scalar, representing the ENM energy at the minimum
#'
#' @return A protein object identical to input, but with v0 added
#' @export
#'
#' @examples
#'
#'@family enm builders
#'
enm_v0_add <- function(prot, v0) {
  graph <- prot$enm$graph
  n_edges <- nrow(graph)
  v0ij <- v0/n_edges
  v0ij <- rep(v0ij,n_edges)
  graph$v0ij <- v0ij
  prot$enm$graph <- graph
  prot
}

#' @rdname enm_v0_add
#' @export
add_v0 <- enm_v0_add
