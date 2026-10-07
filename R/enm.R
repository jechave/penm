
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
  prot <- lst(param = NA, nodes = NA, graph = NA, eij = NA, kmat = NA, nma = NA)
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
  prot$eij <- calculate_enm_edge_geometry(get_xyz(prot), get_graph(prot)$i, get_graph(prot)$j)$eij
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
      mutate(dij = calculate_enm_edge_geometry(xyz, i, j)$dij) %>%
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

#' Calculate vectors, lengths and unit vectors of edges
#'
#' Vectorised over edges: this and `calculate_enm_kmat()` are called at every
#' step of the generalized ENM's minimiser, where looping over edges was
#' measured to make it about four times slower.
#'
#' @param xyz vector of xyz coordinates
#' @param i,j integer vectors of nodes connected in each edge
#' @return a list `(rij, dij, eij)`: `rij` a matrix with n_edges rows and 3
#'   columns (x, y, z), row k the vector from node `i[k]` to node `j[k]`; `dij`
#'   its length; `eij` the unit vector `rij / dij`
#'
#' @family enm builders
#' @noRd
#'
calculate_enm_edge_geometry <- function(xyz, i, j) {
  stopifnot(length(i) == length(j))
  xyz <- my_as_xyz(xyz) # column k is node k
  rij <- t(xyz[, j, drop = FALSE] - xyz[, i, drop = FALSE])
  dij <- sqrt(rowSums(rij^2))
  eij <- rij / dij
  # list(), not lst(): lst() takes ~100 us, half the time of this function,
  # which the minimiser calls at every step
  list(rij = rij, dij = dij, eij = eij)
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


#' Calculate kmat given the ENM graph
#'
#' The Hessian of \eqn{V = \frac12 \sum k_{ij} (d_{ij} - l_{ij})^2} at the
#' conformation whose edge lengths are `graph$dij` and unit vectors `eij`. Each
#' edge contributes the 3 x 3 block
#' \eqn{K_{ij} = -k_{ij} [ e e^T + g_{ij} (I - e e^T) ]}, with
#' \eqn{g_{ij} = (d_{ij} - l_{ij}) / d_{ij}}, at blocks (i, j) and (j, i); the
#' diagonal blocks follow from translational invariance,
#' \eqn{K_{ii} = -\sum_{j \ne i} K_{ij}}.
#'
#' For a network from [set_enm()], `lij` is a copy of `dij`, so `gij` is exactly
#' 0 and each block is \eqn{-k_{ij} e e^T}. An lfenm mutant must never have its
#' kmat recomputed here: in lfenm `K_mut = K_wt` by assumption, while a mutant's
#' `lij != dij` would give it a transverse term.
#'
#' @param graph A tibble representing the ENM graph (with edge information, especially \code{kij}, \code{dij} and \code{lij})
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
#' eij <- calculate_enm_edge_geometry(nodes$xyz, graph$i, graph$j)$eij
#' kmat <- calculate_enm_kmat(graph, eij, nsites = nodes$nsites)
#' }
#'
#' @family enm builders
#' @noRd
#'
#'
calculate_enm_kmat <- function(graph, eij, nsites) {
  stopifnot(max(graph$i, graph$j) <= nsites,
            nrow(graph) == nrow(eij),
            !is.null(graph$dij), !is.null(graph$lij))
  i <- graph$i
  j <- graph$j
  kij <- graph$kij
  gij <- (graph$dij - graph$lij) / graph$dij # relative strain

  # kmat[a, i, b, j] couples coordinate a of node i with coordinate b of node j.
  # Reshaped to 3N x 3N at the end.
  kmat <- array(0, dim = c(3, nsites, 3, nsites))

  # Off-diagonal blocks, element (a, b) of every edge's block at once (vectorised
  # over edges, for the minimiser's sake; see calculate_enm_edge_geometry()). A
  # pair of nodes has at most one edge, so no block is written twice.
  for (a in 1:3) {
    for (b in 1:3) {
      ee_ab <- eij[, a] * eij[, b]
      identity_ab <- as.numeric(a == b)
      kij_ab <- -kij * (ee_ab + gij * (identity_ab - ee_ab))
      kmat[cbind(a, i, b, j)] <- kij_ab
      kmat[cbind(a, j, b, i)] <- kij_ab
    }
  }

  # Diagonal blocks: K_ii = -sum over j of K_ij
  row_sums <- apply(kmat, c(1, 2, 3), sum) # [a, i, b]
  for (site in seq(nsites)) {
    kmat[, site, , site] <- -row_sums[, site, ]
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
#' Given an enm `kmat`, perform NMA: the eigenvalues in ascending order, with
#' the six null modes (rigid-body translations and rotations) left out, and the
#' eigenvectors in the sign convention of [canonical_sign()].
#'
#' The null modes are checked, not assumed: there must be exactly six
#' eigenvalues with \eqn{|\lambda| \le} `null_tol` \eqn{\max |\lambda|}, and
#' none below \eqn{-}`null_tol` \eqn{\max |\lambda|}; otherwise an error. A
#' seventh null eigenvalue means the network is not rigid; a negative one means
#' `kmat` is not at a minimum (e.g. a network with negative springs). The
#' threshold is relative so that it does not depend on the scale of `kij`.
#'
#' @param kmat The K matrix to diagonalize
#' @param null_tol relative threshold for a null eigenvalue
#'
#' @return A list with elements \code{lst(mode,evalue,cmat,umat)}
#'
#' @examples
#' \dontrun{
#' calculate_enm_nma(kmat)
#' }
#'
#'@family enm builders
#' @noRd
#'
#'
calculate_enm_nma <- function(kmat, null_tol = 1e-8) {
  eig <- eigen(kmat, symmetric = TRUE) # eigenvalues in decreasing order
  null_threshold <- null_tol * max(abs(eig$values))

  n_negative <- sum(eig$values < -null_threshold)
  if (n_negative > 0) {
    stop("Hessian has ", n_negative, " negative eigenvalue(s): not a minimum")
  }
  n_null <- sum(abs(eig$values) <= null_threshold)
  if (n_null != 6) {
    stop("Hessian has ", n_null, " null eigenvalues, expected 6: the network is not rigid")
  }

  internal <- which(eig$values > null_threshold)
  ascending <- rev(internal)
  evalue <- eig$values[ascending]
  umat <- canonical_sign(eig$vectors[, ascending, drop = FALSE])
  cmat <- umat %*% ((1 / evalue) * t(umat)) # pseudo-inverse of kmat

  nma <- list(
    mode = seq_along(evalue),
    evalue = evalue,
    cmat = cmat,
    umat = umat
  )
  nma
}
