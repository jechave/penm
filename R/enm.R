
# Create and set prot object ----------------------------------------------


#' Set up 'prot' object
#'
#' @description
#' `set_enm` sets up a `prot` object containing information on ENM structure,
#' parameters, and normal modes. It is the entry point of the package: everything
#' else takes the `prot` it returns.
#'
#' All five arguments are required — there are no defaults.
#'
#' @param pdb   pdb object obtained using [bio3d::read.pdb()]
#' @param node  how network nodes are built: `"ca"` (alpha carbons), `"sc"` (side
#'   chains), or `"cb"` (beta carbons). The long forms `"calpha"`, `"side_chain"`
#'   and `"beta"` are accepted as synonyms.
#' @param model ENM variant, one of `"anm"`, `"ming_wall"`, `"hnm"`, `"hnm0"`,
#'   `"pfanm"`, `"reach"`. These select the spring-constant function applied to each
#'   contact.
#' @param d_max distance cutoff (Å) used to define enm contacts
#' @param frustrated logical value indicating whether to include frustrations in
#'   calculation of kmat. **Only `FALSE` is currently supported**; `TRUE` errors.
#'
#' @returns an object of class `prot`, which is a list `lst(param, nodes, graph, eij, kmat, nma)`
#'
#' @export
#'
#' @seealso [get_mutant_site()] to perturb the result;
#'   [get_prot_property] for the accessors that read a `prot`;
#'   [delta_structure_by_site], [delta_motion_by_site] and [delta_energy] for the
#'   wild-type-vs-mutant comparisons.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, frustrated = FALSE)
#' get_nsites(wt)
#' get_nmodes(wt)
#'
#' # side-chain nodes, a different variant and cutoff
#' wt_sc <- set_enm(pdb_2acy_A, node = "sc", model = "anm",
#'                  d_max = 12.5, frustrated = FALSE)
#' get_nsites(wt_sc)
set_enm <- function(pdb, node, model, d_max, frustrated) {

  stopifnot(!frustrated) # WARNING: need to test frustrated = T option, not sure whether mut_graph is consistent with kmat calculation

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
#' @noRd
#'

create_enm <- function() {
  prot <- lst(param = NA, nodes = NA, graph = NA, eij = NA, kmat = NA, nma = NA)
  class(prot) <- c("prot", class(prot))
  prot
}

#' Set param of prot object
#'
#' @noRd
#'
set_enm_param <- function(prot, node, model, d_max, frustrated) {
  prot$param <- lst(node, model, d_max, frustrated)
  prot
}


#' Set nodes of prot object
#'
#' @noRd
#'
set_enm_nodes <- function(prot, pdb) {
  prot$nodes <- calculate_enm_nodes(pdb, get_enm_node(prot))
  return(prot)
}


#' Set graph of prot object
#'
#' @noRd
#'
set_enm_graph <- function(prot) {
  prot$graph <- calculate_enm_graph(get_xyz(prot), get_pdb_site(prot), get_enm_model(prot), get_d_max(prot))
  prot
}

#' Set eij unit vectors of prot object
#'
#' @noRd
#'
set_enm_eij <- function(prot) {
  prot$eij <- calculate_enm_eij(get_xyz(prot), get_graph(prot)$i, get_graph(prot)$j)
  prot
}

#' Set enm's kmat of prot object
#'
#' @noRd
#'
set_enm_kmat <- function(prot) {
  prot$kmat <- calculate_enm_kmat(get_graph(prot), get_eij(prot), get_nsites(prot), get_frustrated(prot))
  prot
}

#' Set normal-mode-analysis component of prot object
#'
#' @noRd
#'
set_enm_nma <- function(prot) {
  prot$nma <- calculate_enm_nma(get_kmat(prot))
  prot
}


# Calculate prot components -----------------------------------------------


#' Calculate nodes of prot object
#'
#' @param pdb pdb object obtained using bio3d::read.pdb()
#' @param node type, either "ca" or "sc"
#'
#' @returns a list of node properties:  \code{lst(nsites, site, pdb_site, bfactor, xyz)}
#'
#'@family enm builders
#' @noRd
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
  if (node == "cb" | node == "beta") {
    nodes <- prot_cb(pdb)
    return(nodes)
  }
  stop("Error: node must be ca, calpha, sc, side_chain, cb, or beta")
}


#' Calculate ENM graph
#'
#' Calculates graph representation of Elastic Network Model (ENM), the typical relaxed case (lij = dij)
#'
#' @param xyz matrix of size \code{c(3,N)} containing each column the \code{x, y, z} coordinates of each of N nodes
#' @param pdb_site integer vector of size N containing the number of each node (pdb residue number)
#' @param model  character variable specifying the ENM model variant, one of
#'     \code{anm, ming_wall, hnm, hnm0, pfanm, reach}.
#' @param d_max distance-cutoff to define network contacts
#' @return a tibble that contains the graph representation of the network
#'
#' @examples
#' \dontrun{
#'  calculate_enm_graph(xyz, pdb_site, model, d_max)
#' }
#'
#'@family enm builders
#' @noRd
#'
calculate_enm_graph <- function(xyz, pdb_site, model, d_max, ...) {
    # Calculate (relaxed) enm graph from xyz
    # Returns graph for the relaxed case

    # put xyz in the right format and check size
    xyz <- my_as_xyz(xyz)
    nsites <- length(pdb_site)
    stopifnot(ncol(xyz) == nsites)

    # set function to calculate i-j spring constants
    kij_fun <- match.fun(paste0("kij_", model))
    kij_par <- lst(d_max = d_max)

    site <- seq(nsites)
    # calculate graph
    graph <- as_tibble(expand_grid(i = site, j = site)) %>%
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
#' @noRd
#'
dij_edge <- function(xyz, i, j) {
  stopifnot(length(i) == length(j))
  xyz <- my_as_xyz(xyz)
  dij <- rep(NA, length(i))
  for (k in seq(length(i)))  {
    ik <- i[k]
    jk <- j[k]
    rij <- xyz[, jk] - xyz[, ik]
    dij[k] <- sqrt(sum(rij^2))
  }
  dij
}

#' Calculate edge sequence distance
#'
#' @noRd
#'
sdij_edge <- function(pdb_site, i, j) {
  # sequence distance
  stopifnot(length(i) == length(j))
  sdij <- abs(pdb_site[j] - pdb_site[i])
  sdij
}


#' Calculate unit vectors of edges
#'
#' @param i,j integer vectors of nodes connected in each edge
#' @param xyz vector of xyz coordinates
#' @return matrix with n_edge rows and 3 columns (x, y, z)
#'
#' @family enm builders
#' @noRd
#'
calculate_enm_eij <- function(xyz, i, j) {

  stopifnot(length(i) == length(j))
  n_edges <- length(i)
  xyz <- my_as_xyz(xyz)
  # eij <- tibble(eij_x = rep(NA,n_edges),
  #               eij_y = rep(NA,n_edges),
  #               eij_z = rep(NA,n_edges))

  eij <- matrix(NA, n_edges, 3)
  for (k in seq(n_edges)) {
    ik <- i[k]
    jk <- j[k]
    rij <- xyz[, jk] - xyz[, ik]
    eij[k, ] <- rij / sqrt(sum(rij^2))
  }
  eij
}


#' Calculate kmat given the ENM graph
#'
#' @param graph A tibble representing the ENM graph (with edge information, especially \code{kij}
#' @param eij A matrix of size \code{n_edges x 3} of \code{eij} versors directed along ENM contacts
#' @param nsites The number of nodes of the ENM network
#' @param frustrated Logical indicating whether to add frustration or not before calculating \code{kmat}
#'
#' @return The \code{3 nsites x 3 nsites} stiffness matrix of the ENM
#'
#' @examples
#' \dontrun{
#' pdb <- bio3d::read.pdb("2acy")
#' nodes <- calculate_enm_nodes(pdb, node = "ca")
#' graph <- calculate_enm_graph(nodes$xyz, nodes$pdb_site, model = "anm", d_max = 10.5)
#' eij <- calculate_enm_eij(nodes$xyz, graph$i, graph$j)
#' kmat <- calculate_enm_kmat(graph, eij, nsites = nodes$nsites, frustrated = FALSE)
#' }
#'
#' @family enm builders
#' @noRd
#'
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
      gij <- graph$lij[[edge]] / graph$dij[[edge]] - 1
    } else {
      gij <- 0
    }
    eij_v <- eij[edge, ]
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


#' Fix the sign convention of a matrix of eigenvectors
#'
#' `eigen()` determines each eigenvector only up to a factor of -1, and which
#' sign comes back depends on the LAPACK build — so the same `kmat` yields
#' different `umat` on different machines. That makes any stored eigenvector,
#' or any comparison against one, non-portable.
#'
#' The convention adopted here: scale each column so that its largest-magnitude
#' element is positive. This is deterministic given the eigenvector, and leaves
#' every derived quantity unchanged, since all of them (`cmat`, msf, the
#' `delta_*` profiles) are quadratic in `umat`.
#'
#' Note this does not disambiguate degenerate eigenvalues, where any rotation
#' within the degenerate subspace is a valid eigenbasis. ENM spectra of real
#' proteins are generically non-degenerate, so in practice the sign is the whole
#' ambiguity — but a rotation, if one occurred, would survive this.
#'
#' @param umat a matrix whose columns are eigenvectors
#' @return `umat` with each column's largest-magnitude element made positive
#'
#' @noRd
#'
canonical_sign <- function(umat) {
  umat <- as.matrix(umat)
  # sign of the largest-magnitude entry of each column
  pivot <- apply(umat, 2, function(u) u[which.max(abs(u))])
  s <- sign(pivot)
  # a genuinely all-zero column has no sign to fix; leave it be
  s[s == 0] <- 1
  sweep(umat, 2, s, `*`)
}


#' Perform Normal Mode Analysis
#'
#' Given an enm `kmat`, perform NMA
#'
#' @param kmat The K matrix to diagonalize
#' @param too_small=1.e-5 A small value, eigenvectors with eigenvalues larger than `too_small` are discarded
#'
#' @return A list with elements \code{lst(mode,evalue,cmat,umat)}
#'
#' @examples
#' \dontrun{
#' calculate_enm_anm(kmat, too_small = 1.e-10)
#' }
#'
#'@family enm builders
#' @noRd
#'
#'
calculate_enm_nma <- function(kmat, too_small = 1.e-5) {
  eig <- eigen(kmat, symmetric = TRUE)
  evalue <- eig$values
  umat <- eig$vectors
  modes <- evalue > too_small
  evalue <- evalue[modes]
  umat  <- umat[, modes]

  nmodes <- sum(modes)
  mode <- order(seq(nmodes), decreasing = T)
  evalue <- evalue[mode]
  umat <- umat[, mode]
  mode <- mode[mode]

  umat <- canonical_sign(umat)

  cmat <-  umat %*% ((1 / evalue) * t(umat))

  nma <- list(
    mode = mode,
    evalue = evalue,
    cmat = cmat,
    umat = umat
  )
  nma
}
