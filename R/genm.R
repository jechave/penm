# Generalized ENM (genm) -------------------------------------------------------
#
# An ENM whose force constants depend on the equilibrium lengths, kij = k(lij),
# and whose equilibrium lengths depend on the sequence:
#
#   lij = l0ij + delta(i, s_i)ij + delta(j, s_j)ij
#
# l0ij is the edge's length in the pdb, and s_i the allele at site i (allele 0 is
# the pdb's residue, and contributes nothing). The protein is the minimum of
#
#   V(r) = 1/2 sum_ij kij (dij(r) - lij)^2,
#
# with its kmat (the Hessian of V there) and normal modes.
#
# Two objects:
#   - "genm": the parameters (build_enm_from_pdb, genm_mutate)
#   - "genm_prot": the protein they imply (prot_from_enm, genm_add_nma)


# Parameters -------------------------------------------------------------------

#' Build a generalized ENM from a pdb structure
#'
#' One edge for every pair of nodes closer than `d_max_pairs` in the pdb, and
#' one for every i,i+1 pair. `l0ij` is the edge's length in the pdb. Every site
#' starts at allele 0, so `lij = l0ij`, and `kij = k(lij)`.
#'
#' `d_max_pairs` is not the contact cutoff: it must reach out to where `k` is
#' negligible, so that an edge whose `lij` shortens past `d_max` can become a
#' contact. A warning is given when it does not (see
#' [genm_check_d_max_pairs()]).
#'
#' @param pdb pdb object obtained using [bio3d::read.pdb()]
#' @param node `"ca"`, `"sc"` or `"cb"`
#' @param model name of a `kij_*` function, e.g. `"ming_wall"`
#' @param d_max cutoff passed to the `kij_*` function
#' @param d_max_pairs distance (A) in the pdb within which a pair of nodes gets
#'   an edge
#' @param ... further named parameters of the `kij_*` function, e.g. `w`
#' @param ensemble which realization of the mutational process the alleles
#'   refer to; see `?penm_ensemble`
#' @param n_alleles number of alleles per site, including the pdb's (allele 0)
#' @param mut_dl_sigma standard deviation of the change an allele makes to
#'   `lij`
#' @param mut_sd_min edges joining sites less than `mut_sd_min` apart in
#'   sequence are not changed by mutations
#'
#' @returns an object of class `"genm"`: `list(param, nodes, sequence, graph)`,
#'   with `graph` a tibble `(i, j, sdij, l0ij, lij, kij)`
#'
#' @noRd
#'
build_enm_from_pdb <- function(pdb, node, model, d_max, d_max_pairs, ...,
                               ensemble = 1L, n_alleles = 10L,
                               mut_dl_sigma = 0.3, mut_sd_min = 2L) {
  if (d_max_pairs < d_max) stop("d_max_pairs must not be smaller than d_max")
  check_ensemble(ensemble)
  stopifnot(n_alleles >= 2, mut_dl_sigma > 0, mut_sd_min >= 1)

  kij_fun <- match.fun(paste0("kij_", model))
  kij_par <- list(...)
  if (length(kij_par) > 0 && (is.null(names(kij_par)) || any(names(kij_par) == ""))) {
    stop("parameters in ... must be named")
  }
  unknown <- setdiff(names(kij_par), names(formals(kij_fun)))
  if (length(unknown) > 0) {
    stop("kij_", model, " does not take parameter(s): ", paste(unknown, collapse = ", "))
  }

  param <- list(node = node, model = model, d_max = d_max,
                d_max_pairs = d_max_pairs, kij_par = kij_par,
                ensemble = ensemble, n_alleles = as.integer(n_alleles),
                mut_dl_sigma = mut_dl_sigma, mut_sd_min = as.integer(mut_sd_min))

  nodes <- calculate_enm_nodes(pdb, node)
  nsites <- nodes$nsites

  # edges: every pair within d_max_pairs in the pdb, and every i,i+1 pair
  distance <- as.matrix(stats::dist(t(matrix(nodes$xyz, nrow = 3))))
  sequence_distance <- abs(outer(nodes$pdb_site, nodes$pdb_site, "-"))
  is_edge <- upper.tri(distance) & (distance <= d_max_pairs | sequence_distance == 1)
  edges <- which(is_edge, arr.ind = TRUE)
  edges <- edges[order(edges[, 1], edges[, 2]), , drop = FALSE]
  i <- unname(edges[, 1])
  j <- unname(edges[, 2])

  graph <- tibble(
    i = i,
    j = j,
    sdij = sdij_edge(nodes$pdb_site, i, j),
    l0ij = dij_edge(nodes$xyz, i, j)
  )
  graph$lij <- graph$l0ij
  graph$kij <- genm_kij(param, graph$lij, graph$sdij)

  enm <- list(
    param = param,
    nodes = list(nsites = nsites, site = nodes$site,
                 pdb_site = nodes$pdb_site, bfactor = nodes$bfactor),
    sequence = integer(nsites),
    graph = graph
  )
  class(enm) <- c("genm", "list")

  genm_check_d_max_pairs(enm)
  enm
}


#' Mutate a site of a generalized ENM
#'
#' Sets the allele at `site` and recomputes `lij` and `kij`. Mutating back to
#' an earlier allele restores the earlier parameters exactly.
#'
#' @param enm a `"genm"` object
#' @param site the site to mutate (sequential index, not pdb numbering)
#' @param allele the new allele, in `0 .. n_alleles - 1`, different from the
#'   current one
#'
#' @returns the mutant `"genm"` object
#'
#' @noRd
#'
genm_mutate <- function(enm, site, allele) {
  stopifnot(inherits(enm, "genm"))
  if (!(site %in% seq_len(enm$nodes$nsites))) {
    stop("site must be one of 1..", enm$nodes$nsites)
  }
  if (!(allele %in% 0:(enm$param$n_alleles - 1))) {
    stop("allele must be one of 0..", enm$param$n_alleles - 1)
  }
  if (allele == enm$sequence[site]) {
    stop("site ", site, " already has allele ", allele, ": nothing to mutate")
  }

  enm$sequence[site] <- as.integer(allele)
  lij <- genm_lij(enm)
  if (any(lij <= 0)) {
    stop("allele ", allele, " at site ", site, " would make ", sum(lij <= 0),
         " equilibrium length(s) <= 0")
  }
  enm$graph$lij <- lij
  enm$graph$kij <- genm_kij(enm$param, enm$graph$lij, enm$graph$sdij)
  enm
}


#' Equilibrium lengths implied by the sequence
#'
#' \eqn{l_{ij} = l^0_{ij} + \delta(i, s_i)_{ij} + \delta(j, s_j)_{ij}} for every
#' edge, with \eqn{\delta} given by [genm_allele_delta_lij()] for edges with
#' `sdij >= mut_sd_min`, and 0 for the others. Computed from `l0ij` and the
#' sequence only, so equal sequences give identical lengths.
#'
#' @param enm a `"genm"` object
#'
#' @returns `lij`, one per edge
#'
#' @noRd
#'
genm_lij <- function(enm) {
  graph <- enm$graph
  perturbed <- graph$sdij >= enm$param$mut_sd_min
  delta_from_i <- numeric(nrow(graph))
  delta_from_j <- numeric(nrow(graph))

  for (site in which(enm$sequence != 0)) {
    delta_lij <- genm_allele_delta_lij(enm, site, enm$sequence[site])
    site_is_i <- perturbed & graph$i == site
    site_is_j <- perturbed & graph$j == site
    delta_from_i[site_is_i] <- delta_lij[graph$j[site_is_i]]
    delta_from_j[site_is_j] <- delta_lij[graph$i[site_is_j]]
  }

  # always added in this order: floating-point addition is not associative
  (graph$l0ij + delta_from_i) + delta_from_j
}


#' The change an allele makes to the edges of its site
#'
#' A vector of `nsites` normal draws with sd `mut_dl_sigma`, seeded by
#' `(ensemble, site, allele)` without disturbing the caller's RNG. Element `k`
#' is the change to the edge between `site` and site `k` (element `site` itself
#' is not used). Allele 0 changes nothing.
#'
#' @param enm a `"genm"` object
#' @param site sequential site index
#' @param allele an allele in `0 .. n_alleles - 1`
#'
#' @returns a vector of length `nsites`
#'
#' @noRd
#'
genm_allele_delta_lij <- function(enm, site, allele) {
  nsites <- enm$nodes$nsites
  if (allele == 0) return(numeric(nsites))
  seed <- mut_seed(enm$param$ensemble, site, allele)
  with_mut_seed(seed, stats::rnorm(nsites, mean = 0, sd = enm$param$mut_dl_sigma))
}


#' Evaluate kij = k(lij)
#'
#' @param param the `param` list of a `"genm"` object
#' @param lij equilibrium lengths
#' @param sdij sequence separations
#'
#' @returns `kij`; an error if any is not finite or is negative
#'
#' @noRd
#'
genm_kij <- function(param, lij, sdij) {
  kij_fun <- match.fun(paste0("kij_", param$model))
  # the kij_* functions name their first argument dij; here it receives lij
  kij <- do.call(kij_fun, c(list(lij, sdij = sdij, d_max = param$d_max), param$kij_par))
  if (any(!is.finite(kij))) stop("kij_", param$model, " returned non-finite values")
  if (any(kij < 0)) {
    stop("kij_", param$model, " is negative for ", sum(kij < 0), " edge(s): k must be >= 0")
  }
  kij
}


#' Warn when d_max_pairs truncates k
#'
#' Pairs farther apart than `d_max_pairs` in the pdb get no edge, so they can
#' never become contacts. Warns when `k(d_max_pairs)`, for a pair far apart in
#' sequence, is more than `ratio_max` of the largest `k` among edges more than
#' three apart in sequence.
#'
#' @param enm a `"genm"` object
#' @param ratio_max largest acceptable ratio
#'
#' @noRd
#'
genm_check_d_max_pairs <- function(enm, ratio_max = 0.01) {
  param <- enm$param
  far_in_sequence <- 1e6
  k_far <- genm_kij(param, param$d_max_pairs, sdij = far_in_sequence)
  k_near <- max(enm$graph$kij[enm$graph$sdij > 3])
  if (k_far > ratio_max * k_near) {
    warning("k(d_max_pairs = ", param$d_max_pairs, ") is ", signif(k_far / k_near, 2),
            " of the largest non-bonded k for model '", param$model, "'. ",
            "Pairs beyond d_max_pairs get no edge, so the truncation acts as a hard ",
            "cutoff and no contact can form beyond it.", call. = FALSE)
  }
  invisible(enm)
}


# Protein ----------------------------------------------------------------------

#' Build the protein implied by a generalized ENM
#'
#' Minimises V starting from `xyz_seed`, superposes the minimum onto
#' `xyz_seed`, and computes kmat there. Errors if the minimisation does not
#' converge, or ends somewhere that is not a minimum. `nma` is left `NA`: add
#' the modes with [genm_add_nma()].
#'
#' @param enm a `"genm"` object
#' @param xyz_seed starting coordinates, a vector of length `3 * nsites`
#' @param gtol convergence threshold on the largest force component
#' @param max_iter maximum number of steps
#'
#' @returns an object of class `"genm_prot"`:
#'   `list(enm, nodes, v_min, kmat, nma, minimization)`
#'
#' @noRd
#'
prot_from_enm <- function(enm, xyz_seed, gtol = 1e-10, max_iter = 100) {
  stopifnot(inherits(enm, "genm"))
  nsites <- enm$nodes$nsites
  xyz_seed <- as.vector(xyz_seed)
  if (length(xyz_seed) != 3 * nsites) stop("xyz_seed must have length 3 * nsites = ", 3 * nsites)

  # edges with k = 0 contribute nothing to V, its gradient or kmat
  graph <- enm$graph[enm$graph$kij > 0, ]

  xyz <- xyz_seed
  iter <- 0
  repeat {
    force <- -genm_gradient(xyz, graph, nsites)
    if (max(abs(force)) < gtol) break
    if (iter == max_iter) {
      stop("no convergence in ", max_iter, " steps: max |force| = ", signif(max(abs(force)), 3))
    }
    # Newton step: the linear response to the force, dxyz = kmat^-1 force
    kmat <- genm_kmat(xyz, graph, nsites)
    dxyz <- as.vector(solve(genm_kmat_without_rigid_motion(kmat, xyz), force))
    # far from the minimum the full step can overshoot: halve it until V
    # decreases (allowing for rounding in V)
    v_now <- genm_v_xyz(xyz, graph)
    halvings <- 0
    while (genm_v_xyz(xyz + dxyz, graph) > v_now * (1 + 1e-12)) {
      dxyz <- dxyz / 2
      halvings <- halvings + 1
      if (halvings > 50) stop("no step along the linear response decreases V")
    }
    xyz <- xyz + dxyz
    iter <- iter + 1
  }
  if (iter > 0) {
    all_coordinates <- seq_along(xyz)
    xyz <- as.vector(bio3d::fit.xyz(fixed = xyz_seed, mobile = xyz,
                                    fixed.inds = all_coordinates,
                                    mobile.inds = all_coordinates))
  }

  kmat <- genm_kmat(xyz, graph, nsites)
  # a minimum: kmat is positive definite on the internal coordinates
  cholesky <- tryCatch(chol(genm_kmat_without_rigid_motion(kmat, xyz)), error = function(e) NULL)
  if (is.null(cholesky)) stop("the minimisation ended at a point that is not a minimum")

  prot <- list(
    enm = enm,
    nodes = c(enm$nodes, list(xyz = xyz)),
    v_min = genm_v_xyz(xyz, graph),
    kmat = kmat,
    nma = NA,
    minimization = list(iter = iter, force_max = max(abs(force)))
  )
  class(prot) <- c("genm_prot", "list")
  prot
}


#' Add normal modes to a genm_prot
#'
#' @param prot a `"genm_prot"` object
#' @returns `prot`, with `nma` from [calculate_enm_nma()]
#'
#' @noRd
#'
genm_add_nma <- function(prot) {
  stopifnot(inherits(prot, "genm_prot"))
  prot$nma <- calculate_enm_nma(prot$kmat)
  prot
}


#' Superpose a genm_prot onto target coordinates
#'
#' Rotates and translates the protein onto `target`, and recomputes kmat, and
#' the modes if it has them, in the new orientation.
#'
#' @param prot a `"genm_prot"` object
#' @param target coordinates, a vector of length `3 * nsites`
#'
#' @returns `prot`, superposed onto `target`
#'
#' @noRd
#'
genm_superpose_prot <- function(prot, target) {
  stopifnot(inherits(prot, "genm_prot"))
  nsites <- prot$enm$nodes$nsites
  target <- as.vector(target)
  if (length(target) != 3 * nsites) stop("target must have length 3 * nsites = ", 3 * nsites)

  all_coordinates <- seq_along(target)
  xyz <- as.vector(bio3d::fit.xyz(fixed = target, mobile = prot$nodes$xyz,
                                  fixed.inds = all_coordinates,
                                  mobile.inds = all_coordinates))
  graph <- prot$enm$graph[prot$enm$graph$kij > 0, ]
  prot$nodes$xyz <- xyz
  prot$kmat <- genm_kmat(xyz, graph, nsites)
  has_modes <- !identical(prot$nma, NA)
  if (has_modes) prot$nma <- calculate_enm_nma(prot$kmat)
  prot
}


# V, its gradient, and kmat ----------------------------------------------------

#' V at a conformation
#'
#' [v_dij()] with the edge lengths of `xyz`.
#'
#' @param xyz coordinates, a vector of length `3 * nsites`
#' @param graph tibble with columns `i, j, lij, kij`
#'
#' @noRd
#'
genm_v_xyz <- function(xyz, graph) {
  dij <- dij_edge(xyz, graph$i, graph$j)
  v0ij <- 0
  v <- v_dij(dij, v0ij, graph$kij, graph$lij)
  v
}


#' Gradient of V at a conformation
#'
#' \eqn{\partial V / \partial r_i = -\sum_j k_{ij} (d_{ij} - l_{ij}) e_{ij}},
#' with \eqn{e_{ij}} the unit vector from `i` to `j`.
#'
#' @param xyz coordinates, a vector of length `3 * nsites`
#' @param graph tibble with columns `i, j, lij, kij`
#' @param nsites number of nodes
#'
#' @returns a vector of length `3 * nsites`
#'
#' @noRd
#'
genm_gradient <- function(xyz, graph, nsites) {
  dij <- dij_edge(xyz, graph$i, graph$j)
  eij <- calculate_enm_eij(xyz, graph$i, graph$j)
  edge_term <- graph$kij * (dij - graph$lij) * eij # one row per edge

  # each edge adds -edge_term to node i and +edge_term to node j
  n_edges <- nrow(graph)
  incidence <- matrix(0, nrow = n_edges, ncol = nsites)
  incidence[cbind(seq_len(n_edges), graph$i)] <- -1
  incidence[cbind(seq_len(n_edges), graph$j)] <- 1
  gradient <- crossprod(incidence, edge_term) # one row per node

  as.vector(t(gradient))
}


#' kmat at a conformation
#'
#' The Hessian of V. Each edge contributes the 3 x 3 block
#' \eqn{K_{ij} = -k_{ij} [ e e^T + g_{ij} (I - e e^T) ]}, with
#' \eqn{g_{ij} = (d_{ij} - l_{ij}) / d_{ij}}, at blocks (i, j) and (j, i); the
#' diagonal blocks are \eqn{K_{ii} = -\sum_{j \ne i} K_{ij}}. The \eqn{g} term
#' vanishes where `dij = lij`.
#'
#' @param xyz coordinates, a vector of length `3 * nsites`
#' @param graph tibble with columns `i, j, lij, kij`
#' @param nsites number of nodes
#'
#' @returns the `3 nsites x 3 nsites` kmat
#'
#' @noRd
#'
genm_kmat <- function(xyz, graph, nsites) {
  dij <- dij_edge(xyz, graph$i, graph$j)
  eij <- calculate_enm_eij(xyz, graph$i, graph$j)
  gij <- (dij - graph$lij) / dij
  i <- graph$i
  j <- graph$j
  kij <- graph$kij

  # kmat[a, i, b, j] couples coordinate a of node i with coordinate b of node j
  kmat <- array(0, dim = c(3, nsites, 3, nsites))
  for (a in 1:3) {
    for (b in 1:3) {
      ee_ab <- eij[, a] * eij[, b]
      identity_ab <- as.numeric(a == b)
      kij_ab <- -kij * (ee_ab + gij * (identity_ab - ee_ab))
      kmat[cbind(a, i, b, j)] <- kij_ab # element (a, b) of every edge's block
      kmat[cbind(a, j, b, i)] <- kij_ab
    }
  }
  row_sums <- apply(kmat, c(1, 2, 3), sum)
  for (site in seq_len(nsites)) {
    kmat[, site, , site] <- -row_sums[, site, ]
  }

  dim(kmat) <- c(3 * nsites, 3 * nsites)
  kmat
}


#' kmat made stiff against rigid-body motion
#'
#' kmat plus a stiffness, the mean of its diagonal, along each of the six
#' rigid-body motions of `xyz` (three translations, three rotations about the
#' centroid). At a minimum those are kmat's null modes; away from one, the
#' rotations are soft but not null. Either way the result is invertible, and
#' unchanged for any displacement without rigid-body motion.
#'
#' @param kmat a kmat
#' @param xyz the coordinates it was computed at
#'
#' @returns a matrix the size of `kmat`
#'
#' @noRd
#'
genm_kmat_without_rigid_motion <- function(kmat, xyz) {
  position <- matrix(xyz, nrow = 3)            # column k is node k
  r <- position - rowMeans(position)           # relative to the centroid
  nsites <- ncol(r)

  # each motion as a 3 x nsites matrix, one column per node; a rotation about
  # axis u moves the node at r by u x r
  motions <- list(
    move_x = matrix(c(1, 0, 0), nrow = 3, ncol = nsites),
    move_y = matrix(c(0, 1, 0), nrow = 3, ncol = nsites),
    move_z = matrix(c(0, 0, 1), nrow = 3, ncol = nsites),
    turn_x = rbind(0, -r[3, ], r[2, ]),
    turn_y = rbind(r[3, ], 0, -r[1, ]),
    turn_z = rbind(-r[2, ], r[1, ], 0)
  )
  rigid <- qr.Q(qr(sapply(motions, as.vector))) # 3N x 6, orthonormal columns

  stiffness <- mean(diag(kmat))
  kmat + stiffness * tcrossprod(rigid)
}
