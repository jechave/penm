
# Create and set prot object ----------------------------------------------


#' Set up 'prot' object
#'
#' @description
#' `set_enm` sets up a `prot` object containing information on ENM structure,
#' parameters, and normal modes. It is the entry point of the package: everything
#' else takes the `prot` it returns.
#'
#' The protein also carries how it mutates: the mutational model and its
#' parameters are stored with it, and [get_mutant_site()] uses them. A mutant
#' inherits them, so that it can be mutated in turn the same way.
#'
#' `pdb`, `node`, `model` and `d_max` are required. The other arguments have
#' defaults.
#'
#' @details
#' **Sequences and mutations.** Each site carries a letter of an alphabet of
#' `alphabet_size` letters, numbered 0 to `alphabet_size - 1`. In the protein
#' `set_enm()` builds, every site carries letter 0: the residue in the pdb, the
#' one the structure is relaxed for. Its springs have the lengths of the pdb,
#' \eqn{l^0_{ij}}{l0_ij}.
#'
#' Every other letter \eqn{a}{a} at site \eqn{i}{i} shifts the equilibrium
#' length of each of the site's springs, to each site \eqn{j}{j}, by an amount
#' \eqn{\delta_{ij}(a)}{delta_ij(a)}, normal with mean 0 and standard deviation
#' `mut_dl_sigma`. Each amount is drawn once, fixed by `ensemble`, the site and
#' the letter (see `?penm_ensemble`). Letter 0 shifts nothing. A spring's
#' equilibrium length is
#' \deqn{l_{ij} = l^0_{ij} + \delta_{ij}(s_i) + \delta_{ji}(s_j),}{l_ij = l0_ij + delta_ij(s_i) + delta_ji(s_j),}
#' with \eqn{s_i}{s_i} and \eqn{s_j}{s_j} the letters at its two ends. Springs
#' between sites fewer than `mut_sd_min` apart in sequence keep
#' \eqn{l^0_{ij}}{l0_ij}.
#'
#' A mutation changes the letter at one site. From the pdb's letter to another,
#' it changes each of the site's spring lengths by a normal amount with standard
#' deviation `mut_dl_sigma`; between two other letters, by the difference of
#' two such amounts, with standard deviation \eqn{\sqrt{2}}{sqrt(2)} times
#' `mut_dl_sigma`. The
#' lengths depend only on the letters, so mutating a site back restores its
#' earlier lengths. Letters are labels within a site: the same letter at two
#' sites has unrelated effects.
#'
#' @param pdb   pdb object obtained using [bio3d::read.pdb()]
#' @param node  how network nodes are built: `"ca"` (alpha carbons), `"sc"` (side
#'   chains), or `"cb"` (beta carbons). The long forms `"calpha"`, `"side_chain"`
#'   and `"beta"` are accepted as synonyms.
#' @param model ENM variant, one of `"anm"`, `"ming_wall"`, `"hnm"`, `"hnm0"`,
#'   `"pfanm"`, `"reach"`, `"anm_smooth"`, `"ming_wall_smooth"`. These select the
#'   spring-constant function applied to each contact.
#' @param d_max distance cutoff (Å) used to define enm contacts
#' @param ... further parameters of the spring-constant function, by name; for
#'   example `d_max_width`, which `"anm_smooth"` and `"ming_wall_smooth"` require:
#'   the width (Å) of the switch from contact to no contact around `d_max`.
#' @param mut_model the mutational model, `"lfenm"` or `"genm"`. In both, a
#'   mutation changes the letter at one site, which changes the equilibrium
#'   lengths of the site's springs (see Details). In `"lfenm"` the spring
#'   constants do not change, and the mutant's
#'   structure is the linear response to the resulting force. In `"genm"` the
#'   spring constants depend on the equilibrium lengths, and a mutant's
#'   structure is the minimum of its energy.
#' @param d_max_graph distance (Å) in the pdb within which a pair of nodes
#'   gets an edge (pairs adjacent in sequence always get one). For `"lfenm"` it
#'   must equal `d_max`. For `"genm"` it must reach out to where the spring
#'   constant is negligible, so that pairs that mutations bring closer can
#'   become contacts; a warning is given when it does not.
#' @param ensemble an integer naming which realization of the mutational
#'   process the protein's mutants belong to; see `?penm_ensemble`.
#' @param alphabet_size the size of the model's sequence alphabet: the number
#'   of letters each site can carry, numbered 0 to `alphabet_size - 1`. Letter
#'   0 is the residue in the pdb. See Details.
#' @param mut_dl_sigma the standard deviation (Å) of the shift a letter other
#'   than 0 makes to each of its site's spring lengths, relative to the pdb's
#'   length. See Details.
#' @param mut_sd_min springs between sites fewer than `mut_sd_min` apart in
#'   sequence are not changed by mutations.
#'
#' @returns an object of class `prot`, which is a list
#'   `lst(param, nodes, graph, kmat, nma, internal)`. `param` holds the
#'   arguments other than `pdb`, with those given in `...` under `kij_par`.
#'   `nodes` holds the sites, their coordinates and their `sequence` of
#'   letters. `graph` has one row per spring, with its pdb length `l0ij`,
#'   equilibrium length `lij`, spring constant `kij` and current length `dij`.
#'   `kmat` is the network matrix and `nma` its normal modes. `internal` holds
#'   what penm needs for its own computations and is not meant to be read.
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
#'
#' # a smoothed cutoff, and the mutational parameters set explicitly
#' wt_smooth <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall_smooth",
#'                      d_max = 10.5, d_max_width = 1, ensemble = 3,
#'                      mut_dl_sigma = 0.2, mut_sd_min = 3)
#' get_enm_param(wt_smooth)$mut_dl_sigma
set_enm <- function(pdb, node, model, d_max, ..., mut_model = "lfenm",
                    d_max_graph = d_max, ensemble = 1L, alphabet_size = 10L,
                    mut_dl_sigma = 0.3, mut_sd_min = 2L) {

  prot <- create_enm() %>%
    set_enm_param(node = node, model = model, d_max = d_max, kij_par = list(...),
                  mut_model = mut_model, d_max_graph = d_max_graph,
                  ensemble = ensemble, alphabet_size = alphabet_size,
                  mut_dl_sigma = mut_dl_sigma, mut_sd_min = mut_sd_min) %>%
    set_enm_nodes(pdb = pdb)

  if (mut_model == "lfenm") {
    prot <- prot %>%
      set_enm_sequence() %>%
      set_enm_graph() %>%
      set_enm_eij() %>%
      set_enm_kmat() %>%
      set_enm_nma()
  }
  if (mut_model == "genm") {
    prot <- prot %>%
      set_enm_sequence() %>%
      set_genm_graph() %>%
      set_genm_kmat() %>%
      set_enm_nma()
  }

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
#' Checks the parameters, and stores them. See [set_enm()] for their meaning.
#'
#' @noRd
#'
set_enm_param <- function(prot, node, model, d_max, kij_par, mut_model, d_max_graph,
                          ensemble, alphabet_size, mut_dl_sigma, mut_sd_min) {
  if (!(mut_model %in% c("lfenm", "genm"))) {
    stop("mut_model must be \"lfenm\" or \"genm\", not \"", mut_model, "\"")
  }
  if (d_max_graph < d_max) stop("d_max_graph must not be smaller than d_max")
  if (mut_model == "lfenm" && d_max_graph != d_max) {
    stop("for mut_model = \"lfenm\", d_max_graph must equal d_max")
  }
  check_ensemble(ensemble)
  stopifnot(alphabet_size >= 2, mut_dl_sigma > 0, mut_sd_min >= 1)

  # the spring-constant function, and the further parameters it is given
  kij_fun <- match.fun(paste0("kij_", model))
  if (length(kij_par) > 0 && (is.null(names(kij_par)) || any(names(kij_par) == ""))) {
    stop("parameters in ... must be named")
  }
  unknown <- setdiff(names(kij_par), names(formals(kij_fun)))
  if (length(unknown) > 0) {
    stop("kij_", model, " does not take parameter(s): ", paste(unknown, collapse = ", "))
  }

  prot$param <- lst(node, model, d_max, d_max_graph, kij_par, mut_model,
                    ensemble = as.integer(ensemble),
                    alphabet_size = as.integer(alphabet_size), mut_dl_sigma,
                    mut_sd_min = as.integer(mut_sd_min))
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


#' Set the sequence of a prot object: letter 0, the pdb's residue, at every site
#'
#' @noRd
#'
set_enm_sequence <- function(prot) {
  prot$nodes$sequence <- integer(get_nsites(prot))
  prot
}


#' Set graph of prot object
#'
#' @noRd
#'
set_enm_graph <- function(prot) {
  prot$graph <- calculate_enm_graph(get_xyz(prot), get_pdb_site(prot), get_enm_model(prot), get_d_max(prot),
                                    get_enm_param(prot)$kij_par)
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

#' Set the normal modes of a prot object
#'
#' Computes the normal modes of a protein from its network matrix (`kmat`):
#' eigenvalues, eigenvectors and covariance matrix, as the mode accessors
#' ([get_prot_property]) and the measures built on them read them.
#'
#' A protein from [set_enm()] has its modes, and so does an lfenm mutant, which
#' inherits them. A genm mutant from [get_mutant_site()] does not: its modes
#' cost an eigendecomposition, which a long trajectory need not pay at every
#' step. Call `set_enm_nma()` on the proteins whose modes you need.
#'
#' If the protein already has modes, they are recomputed from `kmat`, which
#' gives the same modes.
#'
#' @param prot a `prot`, with or without normal modes
#'
#' @returns `prot`, with its normal modes
#'
#' @export
#'
#' @seealso [get_mutant_site()]; [get_prot_property] for the mode accessors.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
#'               mut_model = "genm", d_max_graph = 14)
#' mut <- get_mutant_site(wt, site_mut = 80, mutation = 3)
#' mut <- set_enm_nma(mut)
#' get_nmodes(mut)
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
#' @param kij_par further named parameters of the spring-constant function
#' @return a tibble that contains the graph representation of the network
#'
#' @examples
#' \dontrun{
#'  calculate_enm_graph(xyz, pdb_site, model, d_max, kij_par = list())
#' }
#'
#'@family enm builders
#' @noRd
#'
calculate_enm_graph <- function(xyz, pdb_site, model, d_max, kij_par) {
    # Calculate (relaxed) enm graph from xyz
    # Returns graph for the relaxed case

    # put xyz in the right format and check size
    xyz <- my_as_xyz(xyz)
    nsites <- length(pdb_site)
    stopifnot(ncol(xyz) == nsites)

    # set function to calculate i-j spring constants
    kij_fun <- match.fun(paste0("kij_", model))
    kij_par <- c(lst(d_max = d_max), kij_par)

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
      mutate(l0ij = dij) %>%
      dplyr::select(edge, i, j, v0ij, sdij, l0ij, lij, kij, dij)

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
