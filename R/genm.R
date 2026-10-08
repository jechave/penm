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
# A genm protein is a prot, with mut_model = "genm" in its param, the sequence
# in its nodes, and l0ij in its graph.


# The wild type ----------------------------------------------------------------
#
# set_enm(..., mut_model = "genm") builds it with these steps. Every site
# starts at allele 0, so lij = l0ij: every edge is at its rest length in the
# pdb structure, which is therefore the minimum, and kmat and the normal modes
# are computed there.

#' Set the graph of a genm prot
#'
#' @noRd
#'
set_genm_graph <- function(prot) {
  prot$graph <- calculate_genm_graph(get_xyz(prot), get_pdb_site(prot), get_enm_param(prot))
  genm_check_d_max_graph(prot)
  prot
}


#' Set the kmat of a genm prot
#'
#' @noRd
#'
set_genm_kmat <- function(prot) {
  graph <- get_graph(prot)
  # edges with k = 0 contribute nothing to kmat
  graph <- graph[graph$kij > 0, ]
  prot$kmat <- genm_kmat(get_xyz(prot), graph, get_nsites(prot))
  prot
}


#' Calculate the graph of a genm prot
#'
#' One edge for every pair of nodes closer than `d_max_graph` in the pdb, and
#' one for every i,i+1 pair. `l0ij` is the edge's length in the pdb. Every site
#' is at allele 0, so `lij = l0ij`, and `kij = k(lij)`.
#'
#' @param xyz coordinates of the pdb structure, a vector of length `3 * nsites`
#' @param pdb_site the pdb numbering of the sites
#' @param param the `param` list of the prot
#'
#' @returns a tibble `(edge, i, j, v0ij, sdij, l0ij, lij, kij, dij)`
#'
#' @noRd
#'
calculate_genm_graph <- function(xyz, pdb_site, param) {
  # edges: every pair within d_max_graph in the pdb, and every i,i+1 pair
  distance <- as.matrix(stats::dist(t(matrix(xyz, nrow = 3))))
  sequence_distance <- abs(outer(pdb_site, pdb_site, "-"))
  is_edge <- upper.tri(distance) & (distance <= param$d_max_graph | sequence_distance == 1)
  edges <- which(is_edge, arr.ind = TRUE)
  edges <- edges[order(edges[, 1], edges[, 2]), , drop = FALSE]
  i <- unname(edges[, 1])
  j <- unname(edges[, 2])

  graph <- tibble(
    edge = paste(i, j, sep = "-"),
    i = i,
    j = j,
    v0ij = 0,
    sdij = sdij_edge(pdb_site, i, j),
    l0ij = dij_edge(xyz, i, j)
  )
  graph$lij <- graph$l0ij
  graph$kij <- genm_kij(param, graph$lij, graph$sdij)
  graph$dij <- graph$l0ij
  graph
}


# Mutations --------------------------------------------------------------------

#' Mutate a site of a genm protein
#'
#' Sets the allele at `site`, recomputes `lij` and `kij`, and finds the new
#' minimum with [genm_minimize()], starting from the protein's structure.
#' Mutating back to an earlier allele restores the earlier parameters exactly.
#'
#' The mutant has no normal modes (`nma` is `NA`): add them with
#' [set_enm_nma()] when they are needed.
#'
#' @param prot a genm `prot`
#' @param site the site to mutate (sequential index, not pdb numbering)
#' @param allele the new allele, in `0 .. n_alleles - 1`, different from the
#'   current one
#'
#' @returns the mutant `prot`
#'
#' @noRd
#'
genm_mutate <- function(prot, site, allele) {
  if (!identical(get_enm_param(prot)$mut_model, "genm")) {
    stop("prot was not built for mut_model = \"genm\"")
  }
  nsites <- get_nsites(prot)
  n_alleles <- get_enm_param(prot)$n_alleles
  if (!(site %in% seq_len(nsites))) {
    stop("site must be one of 1..", nsites)
  }
  if (!(allele %in% 0:(n_alleles - 1))) {
    stop("allele must be one of 0..", n_alleles - 1)
  }
  if (allele == prot$nodes$sequence[site]) {
    stop("site ", site, " already has allele ", allele, ": nothing to mutate")
  }

  prot$nodes$sequence[site] <- as.integer(allele)
  lij <- genm_lij(prot)
  if (any(lij <= 0)) {
    stop("allele ", allele, " at site ", site, " would make ", sum(lij <= 0),
         " equilibrium length(s) <= 0")
  }
  prot$graph$lij <- lij
  prot$graph$kij <- genm_kij(prot$param, prot$graph$lij, prot$graph$sdij)
  genm_minimize(prot)
}


#' Equilibrium lengths implied by the sequence
#'
#' \eqn{l_{ij} = l^0_{ij} + \delta(i, s_i)_{ij} + \delta(j, s_j)_{ij}} for every
#' edge, with \eqn{\delta} given by [genm_allele_delta_lij()] for edges with
#' `sdij >= mut_sd_min`, and 0 for the others. Computed from `l0ij` and the
#' sequence only, so equal sequences give identical lengths.
#'
#' @param prot a genm `prot`
#'
#' @returns `lij`, one per edge
#'
#' @noRd
#'
genm_lij <- function(prot) {
  graph <- prot$graph
  sequence <- prot$nodes$sequence
  perturbed <- graph$sdij >= prot$param$mut_sd_min
  delta_from_i <- numeric(nrow(graph))
  delta_from_j <- numeric(nrow(graph))

  for (site in which(sequence != 0)) {
    delta_lij <- genm_allele_delta_lij(prot, site, sequence[site])
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
#' @param prot a genm `prot`
#' @param site sequential site index
#' @param allele an allele in `0 .. n_alleles - 1`
#'
#' @returns a vector of length `nsites`
#'
#' @noRd
#'
genm_allele_delta_lij <- function(prot, site, allele) {
  nsites <- get_nsites(prot)
  if (allele == 0) return(numeric(nsites))
  seed <- mut_seed(prot$param$ensemble, site, allele)
  with_mut_seed(seed, stats::rnorm(nsites, mean = 0, sd = prot$param$mut_dl_sigma))
}


#' Evaluate kij = k(lij)
#'
#' @param param the `param` list of a genm `prot`
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


#' Warn when d_max_graph truncates k
#'
#' Pairs farther apart than `d_max_graph` in the pdb get no edge, so they can
#' never become contacts. Warns when `k(d_max_graph)`, for a pair far apart in
#' sequence, is more than `ratio_max` of the largest `k` among edges more than
#' three apart in sequence.
#'
#' @param prot a genm `prot`
#' @param ratio_max largest acceptable ratio
#'
#' @noRd
#'
genm_check_d_max_graph <- function(prot, ratio_max = 0.01) {
  param <- prot$param
  far_in_sequence <- 1e6
  k_far <- genm_kij(param, param$d_max_graph, sdij = far_in_sequence)
  k_near <- max(prot$graph$kij[prot$graph$sdij > 3])
  if (k_far > ratio_max * k_near) {
    warning("k(d_max_graph = ", param$d_max_graph, ") is ", signif(k_far / k_near, 2),
            " of the largest non-bonded k for model '", param$model, "'. ",
            "Pairs beyond d_max_graph get no edge, so the truncation acts as a hard ",
            "cutoff and no contact can form beyond it.", call. = FALSE)
  }
  invisible(prot)
}


# The minimum ------------------------------------------------------------------

#' Move a genm protein to the minimum of its V
#'
#' Minimises V starting from the protein's structure, superposes the minimum
#' onto that structure, and computes `dij` and kmat there. Errors if the
#' minimisation does not converge, or ends somewhere that is not a minimum.
#' `nma` is set to `NA`: the modes of the old structure no longer apply.
#'
#' @param prot a genm `prot`, whose `lij` and `kij` may have changed since its
#'   structure was computed
#' @param gtol convergence threshold on the largest force component
#' @param max_iter maximum number of steps
#'
#' @returns `prot` at its minimum
#'
#' @noRd
#'
genm_minimize <- function(prot, gtol = 1e-10, max_iter = 100) {
  nsites <- get_nsites(prot)
  xyz_seed <- get_xyz(prot)

  # edges with k = 0 contribute nothing to V, its gradient or kmat
  graph <- prot$graph[prot$graph$kij > 0, ]

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

  prot$nodes$xyz <- xyz
  prot$graph$dij <- dij_edge(xyz, prot$graph$i, prot$graph$j)
  prot$kmat <- kmat
  prot$nma <- NA
  prot$internal$minimization <- list(iter = iter, force_max = max(abs(force)))
  prot
}


#' Superpose a genm protein onto target coordinates
#'
#' Rotates and translates the protein onto `target`, and recomputes `dij`,
#' kmat, and the modes if it has them, in the new orientation.
#'
#' @param prot a genm `prot`
#' @param target coordinates, a vector of length `3 * nsites`
#'
#' @returns `prot`, superposed onto `target`
#'
#' @noRd
#'
genm_superpose_prot <- function(prot, target) {
  nsites <- get_nsites(prot)
  target <- as.vector(target)
  if (length(target) != 3 * nsites) stop("target must have length 3 * nsites = ", 3 * nsites)

  all_coordinates <- seq_along(target)
  xyz <- as.vector(bio3d::fit.xyz(fixed = target, mobile = get_xyz(prot),
                                  fixed.inds = all_coordinates,
                                  mobile.inds = all_coordinates))
  graph <- prot$graph[prot$graph$kij > 0, ]
  prot$nodes$xyz <- xyz
  prot$graph$dij <- dij_edge(xyz, prot$graph$i, prot$graph$j)
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
