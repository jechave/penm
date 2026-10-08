
# Create and set prot object ----------------------------------------------


#' Set up 'prot' object
#'
#' @description
#' `set_enm` sets up a `prot` object containing information on ENM structure,
#' parameters, and normal modes. It is the entry point of the package: everything
#' else takes the `prot` it returns.
#'
#' All four arguments are required — there are no defaults.
#'
#' @param pdb   pdb object obtained using [bio3d::read.pdb()]
#' @param node  how network nodes are built: `"ca"` (alpha carbons), `"sc"` (side
#'   chains), or `"cb"` (beta carbons). The long forms `"calpha"`, `"side_chain"`
#'   and `"beta"` are accepted as synonyms.
#' @param model ENM variant, one of `"anm"`, `"ming_wall"`, `"hnm"`, `"hnm0"`,
#'   `"pfanm"`, `"reach"`. These select the spring-constant function applied to each
#'   contact.
#' @param d_max distance cutoff (Å) used to define enm contacts
#'
#' @returns an object of class `prot`, which is a list
#'   `lst(param, nodes, graph, kmat, nma, internal)`. `internal` holds what
#'   penm needs for its own computations and is not meant to be read.
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
#'               d_max = 10.5)
#' get_nsites(wt)
#' get_nmodes(wt)
#'
#' # side-chain nodes, a different variant and cutoff
#' wt_sc <- set_enm(pdb_2acy_A, node = "sc", model = "anm",
#'                  d_max = 12.5)
#' get_nsites(wt_sc)
set_enm <- function(pdb, node, model, d_max) {

  prot <- create_enm() %>%
    set_enm_param(node = node, model = model, d_max = d_max) %>%
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
  prot <- lst(param = NA, nodes = NA, graph = NA, kmat = NA, nma = NA, internal = list())
  class(prot) <- c("prot", class(prot))
  prot
}

#' Set param of prot object
#'
#' @noRd
#'
set_enm_param <- function(prot, node, model, d_max) {
  prot$param <- lst(node, model, d_max)
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
  prot$internal$eij <- calculate_enm_eij(get_xyz(prot), get_graph(prot)$i, get_graph(prot)$j)
  prot
}

#' Set enm's kmat of prot object
#'
#' @noRd
#'
set_enm_kmat <- function(prot) {
  prot$kmat <- calculate_enm_kmat(get_graph(prot), get_eij(prot), get_nsites(prot))
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
  xyz <- my_as_xyz(xyz)                 # column k is node k
  # all edges at once, not a loop: genm's minimiser calls this at every step
  rij <- t(xyz[, j] - xyz[, i])         # one row per edge: node i to node j
  dij <- sqrt(rowSums(rij^2))
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
  xyz <- my_as_xyz(xyz)                 # column k is node k
  # all edges at once, not a loop: genm's minimiser calls this at every step
  rij <- t(xyz[, j] - xyz[, i])         # one row per edge: node i to node j
  dij <- sqrt(rowSums(rij^2))
  eij <- rij / dij
  eij
}


#' Calculate kmat given the ENM graph
#'
#' @param graph A tibble representing the ENM graph (with edge information, especially \code{kij}
#' @param eij A matrix of size \code{n_edges x 3} of \code{eij} versors directed along ENM contacts
#' @param nsites The number of nodes of the ENM network
#'
#' @return The \code{3 nsites x 3 nsites} stiffness matrix of the ENM
#'
#' @examples
#' \dontrun{
#' pdb <- bio3d::read.pdb("2acy")
#' nodes <- calculate_enm_nodes(pdb, node = "ca")
#' graph <- calculate_enm_graph(nodes$xyz, nodes$pdb_site, model = "anm", d_max = 10.5)
#' eij <- calculate_enm_eij(nodes$xyz, graph$i, graph$j)
#' kmat <- calculate_enm_kmat(graph, eij, nsites = nodes$nsites)
#' }
#'
#' @family enm builders
#' @noRd
#'
#'
calculate_enm_kmat <- function(graph, eij, nsites) {
  stopifnot(max(graph$i, graph$j) <= nsites,
            nrow(graph) == nrow(eij))
  kmat <- array(0, dim = c(3, nsites, 3, nsites))
  for (edge in seq(nrow(graph))) {
    i <- graph$i[[edge]]
    j <- graph$j[[edge]]
    kij <- graph$kij[[edge]]
    eij_v <- eij[edge, ]
    eij_mat <- tcrossprod(eij_v, eij_v)
    kij_mat <- -kij * eij_mat
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
#' Scales each column by -1 or 1 so that its largest-magnitude element is
#' positive. Every quantity penm derives from `umat` is quadratic in it, so none
#' changes. Degenerate eigenvalues are not disambiguated: there any rotation
#' within the degenerate subspace is a valid eigenbasis, and it survives this.
#'
#' @param umat a matrix whose columns are eigenvectors
#' @return `umat` with each column's largest-magnitude element made positive
#'
#' @noRd
#'
canonical_sign <- function(umat) {
  # eigen() determines each eigenvector only up to a factor of -1, and which
  # sign comes back depends on the LAPACK build, so the same kmat yields
  # different umat on different machines: any stored eigenvector, or any
  # comparison against one, would be non-portable. ENM spectra of real proteins
  # are generically non-degenerate, so in practice the sign is the whole
  # ambiguity.
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
#' Stops if `kmat` has a negative eigenvalue (not a minimum) or does not have
#' exactly six null eigenvalues (not rigid).
#'
#' @param kmat The K matrix to diagonalize
#' @param null_tol=1.e-8 Eigenvalues smaller in magnitude than `null_tol` times the largest are null, and are discarded
#'
#' @return A list with elements \code{lst(mode,evalue,cmat,umat)}
#'
#' @examples
#' \dontrun{
#' calculate_enm_nma(kmat, null_tol = 1.e-10)
#' }
#'
#'@family enm builders
#' @noRd
#'
#'
calculate_enm_nma <- function(kmat, null_tol = 1.e-8) {
  eig <- eigen(kmat, symmetric = TRUE)
  evalue <- eig$values
  umat <- eig$vectors
  too_small <- null_tol * max(abs(evalue))
  if (any(evalue < -too_small)) {
    stop("kmat has ", sum(evalue < -too_small), " negative eigenvalue(s): not a minimum")
  }
  if (sum(abs(evalue) <= too_small) != 6) {
    stop("kmat has ", sum(abs(evalue) <= too_small), " null eigenvalues, expected 6: the network is not rigid")
  }
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
