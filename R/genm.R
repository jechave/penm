# Generalized ENM ------------------------------------------------------------
#
# The state of a protein is a parameter set: a fixed set of springs, one per
# listed pair of nodes, each with an equilibrium length lij and a force constant
# kij = k(lij). Structure, energy and Hessian are derived from it by minimising
#
#   V(r) = 1/2 sum_ij kij (dij(r) - lij)^2.
#
# See dev/explorations/klfenm/klfenm_report.tex for the derivation.
#
# The parameters are a function of the sequence. Each site carries an allele,
# 0 .. n_alleles - 1, with 0 the residue of the pdb. A spring's length is
#
#   lij = l0ij + delta(i, s_i)ij + delta(j, s_j)ij,
#
# with l0ij its length in the pdb, s_i the allele at site i, and delta(site, 0)
# = 0. delta(site, allele) is a stream of draws seeded by hashing (ensemble,
# site, allele), as for get_mutant_site(); see ?penm_ensemble. The spring from
# site to site k takes the k-th draw, so delta(i, s_i)ij depends on the spring's
# other end, not on which other springs exist. Because lij is always computed
# from l0ij and the sequence, never accumulated, the same sequence gives the
# same parameters bit for bit, however it was reached.
#
# Two objects:
#   - "genm": the parameters (build_enm_from_pdb, genm_mutate)
#   - "genm_prot": what the parameters imply (prot_from_enm; genm_add_nma adds
#     the normal modes on request)
#
# Coordinates are a flat vector of length 3 * nsites, x1 y1 z1 x2 y2 z2 ..., as
# in the rest of penm.
#
# Isolated from the set_enm() / get_mutant_site() machinery on purpose: nothing
# here calls it. Shared with it: calculate_enm_nodes() (pdb parsing), the kij_*
# functions, canonical_sign(), and the mutant-key hashing in R/seed.R.


# Build the parameters ---------------------------------------------------------

#' Build a generalized ENM from a pdb structure
#'
#' The only place where a structure determines parameters. There is one spring
#' for every pair of nodes closer than `d_max_pairs` in the pdb, plus one for
#' every i,i+1 pair regardless of distance; each spring's `l0ij` is its distance
#' in the pdb. Every site starts at allele 0, so `lij = l0ij` and
#' `kij = k(lij)`.
#'
#' `d_max_pairs` is not the contact cutoff. The springs must reach out to where
#' `k` is negligible, so that a spring whose `lij` shortens past `d_max` can
#' become a contact later. A warning is given when `k(d_max_pairs)` is not small
#' compared with `k` at contact distances (see [genm_check_d_max_pairs()]).
#'
#' The last four arguments define the mutational process, and are fixed for
#' every protein derived from this one: they decide what each allele means.
#'
#' @param pdb pdb object obtained using [bio3d::read.pdb()]
#' @param node `"ca"`, `"sc"` or `"cb"` (or `"calpha"`, `"side_chain"`, `"beta"`)
#' @param model name of a `kij_*` function, e.g. `"ming_wall"`, `"anm"`,
#'   `"ming_wall_smooth"`
#' @param d_max cutoff passed to the `kij_*` function
#' @param d_max_pairs distance (A) in the pdb within which a pair of nodes gets a
#'   spring
#' @param ... further named parameters of the `kij_*` function, e.g. `w` for the
#'   smoothed models. Names it does not accept are an error.
#' @param ensemble which realization of the mutational process the alleles refer
#'   to; see `?penm_ensemble`
#' @param n_alleles number of alleles per site, including the pdb's (allele 0)
#' @param mut_dl_sigma standard deviation of the change in `lij` that an allele
#'   contributes to each of its site's springs
#' @param mut_sd_min springs joining sites less than `mut_sd_min` apart in
#'   sequence are not changed by mutations; the default 2 leaves the i,i+1 bonds
#'   alone, 1 perturbs every spring
#'
#' @returns an object of class `"genm"`: `list(param, nodes, sequence, springs)`,
#'   with `sequence` all zeros and `springs` a tibble
#'   `(i, j, sdij, l0ij, lij, kij)` sorted by `(i, j)`
#'
#' @noRd
#'
build_enm_from_pdb <- function(pdb, node, model, d_max, d_max_pairs, ...,
                               ensemble = 1L, n_alleles = 10L,
                               mut_dl_sigma = 0.3, mut_sd_min = 2L) {
  stopifnot(is.numeric(d_max), length(d_max) == 1, d_max > 0,
            is.numeric(d_max_pairs), length(d_max_pairs) == 1)
  if (d_max_pairs < d_max) {
    stop("d_max_pairs (", d_max_pairs, ") must not be smaller than d_max (",
         d_max, "): pairs inside the cutoff would be left out")
  }
  check_ensemble(ensemble)
  genm_check_integer(n_alleles, "n_alleles", min = 2)
  genm_check_integer(mut_sd_min, "mut_sd_min", min = 1)
  if (!is.numeric(mut_dl_sigma) || length(mut_dl_sigma) != 1 ||
      !is.finite(mut_dl_sigma) || mut_dl_sigma <= 0) {
    stop("mut_dl_sigma must be a single positive number")
  }

  param <- list(node = node, model = model, d_max = d_max,
                d_max_pairs = d_max_pairs, kij_par = list(...),
                ensemble = ensemble, n_alleles = as.integer(n_alleles),
                mut_dl_sigma = mut_dl_sigma, mut_sd_min = as.integer(mut_sd_min))
  genm_kij_fun(param) # validate model and kij_par before any work

  nodes <- calculate_enm_nodes(pdb, node)
  nsites <- nodes$nsites

  # Which pairs get a spring: within d_max_pairs in the pdb, or bonded (i,i+1).
  distance <- as.matrix(stats::dist(t(matrix(nodes$xyz, nrow = 3))))
  seq_distance <- abs(outer(nodes$pdb_site, nodes$pdb_site, "-"))
  has_spring <- upper.tri(distance) & (distance <= d_max_pairs | seq_distance == 1)
  pairs <- which(has_spring, arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
  i <- unname(pairs[, 1])
  j <- unname(pairs[, 2])

  springs <- tibble(
    i = i,
    j = j,
    sdij = abs(nodes$pdb_site[j] - nodes$pdb_site[i]),
    # recomputed with genm_geometry(), the distance V uses, so that the wild
    # type sits at V = 0 exactly rather than to within rounding
    l0ij = genm_geometry(nodes$xyz, i, j)$dij
  )
  springs$lij <- springs$l0ij
  springs$kij <- genm_kij(param, springs$lij, springs$sdij)

  enm <- list(
    param = param,
    nodes = list(nsites = nsites, site = nodes$site,
                 pdb_site = nodes$pdb_site, bfactor = nodes$bfactor),
    sequence = integer(nsites),
    springs = springs
  )
  class(enm) <- c("genm", "list")

  genm_check_d_max_pairs(enm)
  enm
}


# Mutate the parameters --------------------------------------------------------

#' Mutate a site of a generalized ENM
#'
#' Changes the allele at `site` from its current value to `allele`, and
#' recomputes `lij` and `kij` of the springs connected to `site`. Produces no
#' structure and no energy, and looks at no other protein.
#'
#' `allele` must differ from the current one: a mutation that changes nothing is
#' an error. Mutating back to an earlier allele, 0 included, is a mutation like
#' any other, and restores that site's springs exactly.
#'
#' @param enm a `"genm"` object
#' @param site the site to mutate (sequential index, not pdb numbering)
#' @param allele the new allele, in `0 .. n_alleles - 1`
#'
#' @returns the mutant `"genm"` object
#'
#' @noRd
#'
genm_mutate <- function(enm, site, allele) {
  stopifnot(inherits(enm, "genm"))
  genm_check_integer(site, "site", min = 1, max = enm$nodes$nsites)
  genm_check_integer(allele, "allele", min = 0, max = enm$param$n_alleles - 1)
  if (allele == enm$sequence[site]) {
    stop("site ", site, " already has allele ", allele, ": nothing to mutate")
  }

  enm$sequence[site] <- as.integer(allele)

  rows <- genm_site_rows(enm, site)
  new_lij <- genm_lij(enm, rows)
  if (any(new_lij <= 0)) {
    stop("allele ", allele, " at site ", site, " would make ", sum(new_lij <= 0),
         " equilibrium length(s) <= 0")
  }
  enm$springs$lij[rows] <- new_lij
  enm$springs$kij <- genm_kij(enm$param, enm$springs$lij, enm$springs$sdij)
  enm
}


#' Equilibrium lengths implied by the sequence
#'
#' For each spring in `rows`, \eqn{l_{ij} = l^0_{ij} + \delta(i, s_i)_{ij} +
#' \delta(j, s_j)_{ij}}: its pdb length, plus what the allele at each of its two
#' ends contributes. Uses `l0ij` and `enm$sequence` only, never the current
#' `lij`, and always adds in this order, so equal sequences give bitwise-equal
#' lengths.
#'
#' @param enm a `"genm"` object
#' @param rows row indices of `enm$springs`
#'
#' @returns the lengths, one per row
#'
#' @noRd
#'
genm_lij <- function(enm, rows) {
  springs <- enm$springs
  dl_from_i <- genm_dl_from_end(enm, rows, springs$i[rows])
  dl_from_j <- genm_dl_from_end(enm, rows, springs$j[rows])
  (springs$l0ij[rows] + dl_from_i) + dl_from_j
}


#' What the alleles at one end of each spring contribute to its length
#'
#' @param enm a `"genm"` object
#' @param rows row indices of `enm$springs`
#' @param end_site for each of `rows`, the site at the end being considered
#'   (`springs$i[rows]` or `springs$j[rows]`)
#'
#' @returns one value per row: \eqn{\delta(s, \text{allele}_s)} for the spring,
#'   with `s = end_site`
#'
#' @noRd
#'
genm_dl_from_end <- function(enm, rows, end_site) {
  dl <- numeric(length(rows))
  mutated <- unique(end_site[enm$sequence[end_site] != 0]) # allele 0 contributes 0
  for (site in mutated) {
    at_this_site <- end_site == site
    site_dl <- genm_site_dl(enm, site, enm$sequence[site])
    position <- match(rows[at_this_site], genm_site_rows(enm, site))
    dl[at_this_site] <- site_dl[position]
  }
  dl
}


#' The change an allele makes to its site's springs
#'
#' One value per spring connected to `site`, in the order of [genm_site_rows()]:
#' 0 for allele 0 and for springs with `sdij < mut_sd_min`; otherwise a normal
#' draw with sd `mut_dl_sigma`, seeded by hashing `(ensemble, site, allele)`
#' without disturbing the caller's RNG.
#'
#' The spring joining `site` to site `k` takes the `k`-th value of that stream.
#' Its change therefore depends on `(ensemble, site, allele, k)` only, not on
#' which other springs exist: networks built with a different `d_max_pairs` or
#' `mut_sd_min` give the springs they share the same change.
#'
#' @param enm a `"genm"` object
#' @param site sequential site index
#' @param allele an allele in `0 .. n_alleles - 1`
#'
#' @noRd
#'
genm_site_dl <- function(enm, site, allele) {
  param <- enm$param
  springs <- enm$springs
  rows <- genm_site_rows(enm, site)
  dl <- numeric(length(rows))
  if (allele == 0) return(dl)

  partner <- ifelse(springs$i[rows] == site, springs$j[rows], springs$i[rows])
  # One draw per site, nsites in all, so that the draw for the spring to site k
  # is simply draws[k]. Not nsites - 1: the value at k = site is never used (no
  # spring joins a site to itself), and skipping it would mean shifting every
  # index above site down by one. Each value depends only on its position in
  # the stream, not on the stream's length.
  draws <- with_mut_seed(
    mut_seed(param$ensemble, site, allele),
    stats::rnorm(enm$nodes$nsites, mean = 0, sd = param$mut_dl_sigma)
  )
  perturbed <- springs$sdij[rows] >= param$mut_sd_min
  dl[perturbed] <- draws[partner[perturbed]]
  dl
}


#' Rows of `enm$springs` connected to a site
#'
#' @param enm a `"genm"` object
#' @param site sequential site index
#'
#' @returns integer row indices, in increasing order
#'
#' @noRd
#'
genm_site_rows <- function(enm, site) {
  which(enm$springs$i == site | enm$springs$j == site)
}


#' Springs connected to a site
#'
#' The springs a mutation at `site` can change (those with
#' `sdij >= mut_sd_min`), together with any it leaves alone.
#'
#' @param enm a `"genm"` object
#' @param site sequential site index
#'
#' @returns a tibble `(pair, i, j, sdij, l0ij, lij, kij)`, `pair` being the row
#'   index in `enm$springs`
#'
#' @noRd
#'
genm_site_springs <- function(enm, site) {
  stopifnot(inherits(enm, "genm"))
  genm_check_integer(site, "site", min = 1, max = enm$nodes$nsites)
  rows <- genm_site_rows(enm, site)
  springs <- enm$springs[rows, ]
  tibble(pair = rows, i = springs$i, j = springs$j, sdij = springs$sdij,
         l0ij = springs$l0ij, lij = springs$lij, kij = springs$kij)
}


# From parameters to protein ---------------------------------------------------

#' Build the protein implied by a generalized ENM
#'
#' Minimises the potential starting from `xyz_seed`, superposes the minimum onto
#' `xyz_seed`, and evaluates the Hessian there, including the transverse term of
#' frustrated springs. `kij` is read from `enm`, never evaluated here.
#'
#' The result depends on `enm` alone when the potential has a single minimum.
#' A frustrated network can have several, and then the one returned is the one
#' reached from `xyz_seed`. `xyz_seed` also fixes the frame: the minimum is
#' superposed onto it.
#'
#' Fails if the minimiser does not converge, or if the stationary point reached
#' is not a minimum of a rigid network.
#'
#' The Hessian is not diagonalised: `nma` is left `NA`, so every getter that
#' reads the modes errors until [genm_add_nma()] fills it. A trajectory of
#' mutations needs only the minimum, not the modes.
#'
#' @param enm a `"genm"` object
#' @param xyz_seed starting coordinates, a vector of length `3 * nsites`
#' @param gtol convergence threshold on the largest gradient component
#' @param max_iter maximum number of Newton steps
#'
#' @returns an object of class `"genm_prot"`:
#'   `list(enm, nodes, v_min, kmat, nma, minimization)`, with `nma = NA`.
#'   `nodes`, `kmat` and (once filled) `nma` have the layout of a `prot` from
#'   [set_enm()], so the getters and analysis functions that read only those
#'   work on it. Those that read `graph` or `param` do not.
#'
#' @noRd
#'
prot_from_enm <- function(enm, xyz_seed, gtol = 1e-10, max_iter = 100) {
  stopifnot(inherits(enm, "genm"))
  nsites <- enm$nodes$nsites
  xyz_seed <- as.vector(xyz_seed)
  if (!is.numeric(xyz_seed) || length(xyz_seed) != 3 * nsites || any(!is.finite(xyz_seed))) {
    stop("xyz_seed must be a finite numeric vector of length 3 * nsites = ", 3 * nsites)
  }

  springs <- genm_springs_with_k(enm)
  minimum <- genm_minimize(springs, nsites, xyz_seed, gtol, max_iter)

  prot <- list(
    enm = enm,
    nodes = c(enm$nodes, list(xyz = minimum$xyz)),
    v_min = genm_energy(minimum$xyz, springs),
    kmat = minimum$kmat,
    nma = NA,
    minimization = list(iter = minimum$iter, grad_max = minimum$grad_max)
  )
  class(prot) <- c("genm_prot", "list")
  prot
}


#' Add normal modes to a genm_prot
#'
#' Diagonalises the Hessian the protein already holds and returns the protein
#' with `nma` filled in; every other component is left as it is. Calling it on a
#' protein that already has modes recomputes them, with the same result: they
#' are a function of `kmat` alone.
#'
#' @param prot a `"genm_prot"` object
#' @returns `prot`, with `nma = list(mode, evalue, cmat, umat)` (see
#'   [genm_nma()])
#'
#' @noRd
#'
genm_add_nma <- function(prot) {
  stopifnot(inherits(prot, "genm_prot"))
  prot$nma <- genm_nma(prot$kmat)
  prot
}


#' Superpose a genm_prot onto target coordinates
#'
#' Rotates and translates the protein so that its structure has the smallest
#' RMSD to `target`, and recomputes everything that depends on orientation from
#' the new coordinates: the Hessian, and the normal modes if the protein has
#' them. The parameters (`enm`) and the minimum energy do not depend on
#' orientation and are kept.
#'
#' Two proteins must be in the same orientation before they are compared site
#' by site or mode by mode; this puts one in the orientation of the other.
#'
#' @param prot a `"genm_prot"` object
#' @param target coordinates to superpose onto, a vector of length
#'   `3 * nsites`, e.g. `get_xyz()` of another protein
#'
#' @returns `prot`, superposed onto `target`
#'
#' @noRd
#'
genm_superpose_prot <- function(prot, target) {
  stopifnot(inherits(prot, "genm_prot"))
  nsites <- prot$enm$nodes$nsites
  target <- as.vector(target)
  if (!is.numeric(target) || length(target) != 3 * nsites || any(!is.finite(target))) {
    stop("target must be a finite numeric vector of length 3 * nsites = ", 3 * nsites)
  }

  xyz <- genm_superpose(prot$nodes$xyz, target)
  prot$nodes$xyz <- xyz
  prot$kmat <- genm_hessian(xyz, genm_springs_with_k(prot$enm), nsites)
  has_modes <- !identical(prot$nma, NA)
  if (has_modes) prot$nma <- genm_nma(prot$kmat)
  prot
}


#' Springs with a nonzero force constant
#'
#' Springs with `k = 0` contribute exactly nothing to V, its gradient or its
#' Hessian; leaving them out is for speed only (about a third of the time for
#' 2acy).
#'
#' @param enm a `"genm"` object
#' @returns the rows of `enm$springs` with `kij > 0`
#'
#' @noRd
#'
genm_springs_with_k <- function(enm) {
  enm$springs[enm$springs$kij > 0, ]
}


#' Minimum energy of a genm_prot
#'
#' @param prot a `"genm_prot"` object
#' @returns \eqn{V(r_e)}, the potential at the minimum
#'
#' @noRd
#'
genm_v_min <- function(prot) {
  stopifnot(inherits(prot, "genm_prot"))
  prot$v_min
}


# Potential, gradient, Hessian -------------------------------------------------

#' Spring vectors and lengths at a conformation
#'
#' @param xyz coordinates, a vector of length `3 * nsites`
#' @param i,j the two nodes of each spring
#'
#' @returns `list(rij, dij)`: `rij` a matrix with one row per spring, the vector
#'   from node `i` to node `j`; `dij` its length
#'
#' @noRd
#'
genm_geometry <- function(xyz, i, j) {
  position <- matrix(xyz, nrow = 3) # column k is node k
  rij <- t(position[, j, drop = FALSE] - position[, i, drop = FALSE])
  list(rij = rij, dij = sqrt(rowSums(rij^2)))
}


#' Potential energy
#'
#' \eqn{V = \frac12 \sum_{ij} k_{ij} (d_{ij} - l_{ij})^2}
#'
#' @param xyz coordinates, a vector of length `3 * nsites`
#' @param springs tibble with columns `i, j, lij, kij`: any set of springs
#'
#' @noRd
#'
genm_energy <- function(xyz, springs) {
  dij <- genm_geometry(xyz, springs$i, springs$j)$dij
  0.5 * sum(springs$kij * (dij - springs$lij)^2)
}


#' Gradient of the potential
#'
#' \eqn{\partial V / \partial r_i = -\sum_j k_{ij} (d_{ij} - l_{ij}) e_{ij}}, with
#' \eqn{e_{ij}} the unit vector from `i` to `j`. Each spring contributes
#' \eqn{-k (d - l) e} to node `i` and \eqn{+k (d - l) e} to node `j`.
#'
#' @param xyz coordinates, a vector of length `3 * nsites`
#' @param springs tibble with columns `i, j, lij, kij`: any set of springs
#' @param nsites number of nodes
#'
#' @returns a vector of length `3 * nsites`
#'
#' @noRd
#'
genm_gradient <- function(xyz, springs, nsites) {
  geometry <- genm_geometry(xyz, springs$i, springs$j)
  stretch <- geometry$dij - springs$lij
  unit <- geometry$rij / geometry$dij
  # k (d - l) e for each spring, one row per spring
  spring_term <- springs$kij * stretch * unit

  gradient <- matrix(0, nrow = 3, ncol = nsites) # column k is node k
  for (a in 1:3) {
    gradient[a, ] <- genm_sum_by_node(-spring_term[, a], springs$i, nsites) +
      genm_sum_by_node(spring_term[, a], springs$j, nsites)
  }
  as.vector(gradient)
}


#' Sum values per node
#'
#' @param values one value per spring
#' @param node the node each value belongs to
#' @param nsites number of nodes
#'
#' @returns a vector of length `nsites`: for each node, the sum of its values
#'   (0 for a node with none)
#'
#' @noRd
#'
genm_sum_by_node <- function(values, node, nsites) {
  total <- numeric(nsites)
  per_node <- rowsum(values, node) # one row per node that has values, named by node
  total[as.integer(rownames(per_node))] <- per_node
  total
}


#' Hessian of the potential
#'
#' Each spring (i, j) contributes the 3 x 3 block
#' \eqn{K_{ij} = -k_{ij} [ e e^T + g_{ij} (I - e e^T) ]}, with
#' \eqn{g_{ij} = (d_{ij} - l_{ij}) / d_{ij}}, at blocks (i, j) and (j, i). The
#' diagonal blocks follow from translational invariance,
#' \eqn{K_{ii} = -\sum_{j \ne i} K_{ij}}. Valid at any conformation, not only at
#' a minimum.
#'
#' @param xyz coordinates, a vector of length `3 * nsites`
#' @param springs tibble with columns `i, j, lij, kij`: any set of springs
#' @param nsites number of nodes
#'
#' @returns the `3 nsites x 3 nsites` Hessian
#'
#' @noRd
#'
genm_hessian <- function(xyz, springs, nsites) {
  geometry <- genm_geometry(xyz, springs$i, springs$j)
  unit <- geometry$rij / geometry$dij              # e, one row per spring
  g <- (geometry$dij - springs$lij) / geometry$dij  # relative strain
  k <- springs$kij

  # The Hessian as a 4-index array: kmat[a, i, b, j] couples coordinate a of
  # node i with coordinate b of node j. Reshaped to 3N x 3N at the end.
  kmat <- array(0, dim = c(3, nsites, 3, nsites))

  # Off-diagonal blocks, element (a, b) for every spring at once. A pair of
  # nodes has at most one spring, so no block is written twice.
  for (a in 1:3) {
    for (b in 1:3) {
      ee_ab <- unit[, a] * unit[, b]
      identity_ab <- as.numeric(a == b)
      block_ab <- -k * (ee_ab + g * (identity_ab - ee_ab))
      kmat[cbind(a, springs$i, b, springs$j)] <- block_ab
      kmat[cbind(a, springs$j, b, springs$i)] <- block_ab
    }
  }

  # Diagonal blocks: K_ii = -sum over j of K_ij
  row_sums <- apply(kmat, c(1, 2, 3), sum) # [a, i, b]
  for (i in seq_len(nsites)) {
    kmat[, i, , i] <- -row_sums[, i, ]
  }

  dim(kmat) <- c(3 * nsites, 3 * nsites)
  kmat
}


# Minimisation ----------------------------------------------------------------

#' Minimise the potential by Newton's method
#'
#' Repeats Newton steps ([genm_newton_step()]) with a backtracking line search
#' ([genm_line_search()]) until the largest gradient component is below `gtol`.
#' Then superposes the result onto `xyz_seed` (unless no step was taken, in
#' which case it is `xyz_seed` itself) and checks that it is a minimum of a
#' rigid network ([genm_is_minimum()]). Errors if either fails.
#'
#' @param springs tibble with columns `i, j, lij, kij`: any set of springs
#' @param nsites number of nodes
#' @param xyz_seed starting coordinates
#' @param gtol convergence threshold on `max(abs(gradient))`
#' @param max_iter maximum number of Newton steps
#'
#' @returns `list(xyz, kmat, iter, grad_max)`
#'
#' @noRd
#'
genm_minimize <- function(springs, nsites, xyz_seed, gtol, max_iter) {
  stopifnot(is.numeric(gtol), gtol > 0, is.numeric(max_iter), max_iter >= 0)
  xyz <- xyz_seed
  iter <- 0
  repeat {
    gradient <- genm_gradient(xyz, springs, nsites)
    grad_max <- max(abs(gradient))
    if (grad_max < gtol) break
    if (iter >= max_iter) {
      stop("genm_minimize did not converge in ", max_iter,
           " iterations: max |gradient| = ", signif(grad_max, 3))
    }
    step <- genm_newton_step(xyz, gradient, springs, nsites)
    xyz <- genm_line_search(xyz, step, gradient, springs, nsites)
    iter <- iter + 1
  }

  # with no step taken xyz is xyz_seed itself; superposing would only add rounding
  if (iter > 0) xyz <- genm_superpose(xyz, xyz_seed)

  kmat <- genm_hessian(xyz, springs, nsites)
  if (!genm_is_minimum(xyz, kmat)) {
    stop("genm_minimize reached a stationary point that is not a minimum of a rigid network ",
         "(the Hessian is not positive definite on the internal coordinates)")
  }
  list(xyz = xyz, kmat = kmat, iter = iter, grad_max = grad_max)
}


#' One Newton step, without rigid-body motion
#'
#' The Newton step solves \eqn{K \, step = -gradient}, but \eqn{K} is singular
#' along the six rigid-body directions (three translations, three rotations).
#' So \eqn{K} is made invertible by adding a stiffness \eqn{c} along exactly
#' those directions, \eqn{K + c P} with \eqn{P} the projector onto them and
#' \eqn{c} the mean of \eqn{K}'s diagonal; the rigid-body part of the resulting
#' step is then discarded.
#'
#' Far from a minimum, compressed springs can make \eqn{K + c P} indefinite, and
#' the plain Newton step may then point uphill. In that case every eigenvalue is
#' replaced by its absolute value, which gives a step that still goes downhill.
#'
#' @param xyz coordinates, a vector of length `3 * nsites`
#' @param gradient the gradient at `xyz`
#' @param springs tibble with columns `i, j, lij, kij`
#' @param nsites number of nodes
#'
#' @returns the step, a vector of length `3 * nsites`
#'
#' @noRd
#'
genm_newton_step <- function(xyz, gradient, springs, nsites) {
  kmat <- genm_hessian(xyz, springs, nsites)
  rigid <- genm_rigid_basis(xyz) # 3N x 6, orthonormal columns
  kmat_shifted <- kmat + mean(diag(kmat)) * tcrossprod(rigid)

  # Cholesky fails exactly when kmat_shifted is not positive definite
  cholesky <- tryCatch(chol(kmat_shifted), error = function(e) NULL)
  if (!is.null(cholesky)) {
    # kmat_shifted = t(cholesky) %*% cholesky; solve in two triangular steps
    step <- -backsolve(cholesky, backsolve(cholesky, gradient, transpose = TRUE))
  } else {
    eig <- eigen(kmat_shifted, symmetric = TRUE)
    # keep tiny eigenvalues away from zero so the step stays finite
    floor <- 1e-8 * max(abs(eig$values))
    stiffness <- pmax(abs(eig$values), floor)
    step <- -eig$vectors %*% (crossprod(eig$vectors, gradient) / stiffness)
  }

  rigid_part <- rigid %*% crossprod(rigid, step)
  as.vector(step - rigid_part)
}


#' Backtracking line search along a Newton step
#'
#' Tries the full step, then half, a quarter, ..., and accepts the first that
#' lowers V enough (the Armijo condition: by at least a small fraction of the
#' decrease the gradient predicts).
#'
#' Close to the minimum, the decrease in V becomes smaller than the rounding
#' error of V itself, and the Armijo condition can no longer be checked. There a
#' step is accepted if V has not risen by more than rounding and the gradient
#' has become smaller.
#'
#' @param xyz current coordinates
#' @param step the proposed step
#' @param gradient the gradient at `xyz`
#' @param springs tibble with columns `i, j, lij, kij`
#' @param nsites number of nodes
#'
#' @returns the new coordinates
#'
#' @noRd
#'
genm_line_search <- function(xyz, step, gradient, springs, nsites) {
  sufficient_decrease <- 1e-4                       # Armijo constant
  smallest_fraction <- 1e-10                        # give up below this
  v_now <- genm_energy(xyz, springs)
  v_rounding <- 64 * .Machine$double.eps * abs(v_now) # rounding error of V
  grad_max_now <- max(abs(gradient))

  slope <- sum(gradient * step) # dV/dt along the step, at t = 0
  if (slope >= 0) stop("genm_minimize: Newton step is not a descent direction")

  fraction <- 1
  repeat {
    xyz_new <- xyz + fraction * step
    v_new <- genm_energy(xyz_new, springs)

    if (v_new <= v_now + sufficient_decrease * fraction * slope) return(xyz_new)

    v_unchanged <- v_new <= v_now + v_rounding
    if (v_unchanged) {
      grad_max_new <- max(abs(genm_gradient(xyz_new, springs, nsites)))
      if (grad_max_new < grad_max_now) return(xyz_new)
    }

    fraction <- fraction / 2
    if (fraction < smallest_fraction) {
      stop("genm_minimize: line search failed to decrease the energy")
    }
  }
}


#' Is a stationary point a minimum of a rigid network?
#'
#' True when the Hessian is positive definite on the internal coordinates,
#' i.e. when \eqn{K + c P} (see [genm_newton_step()]) has a Cholesky
#' factorisation. False for a saddle point, and for a network that is not rigid.
#'
#' @param xyz coordinates of the stationary point
#' @param kmat the Hessian there
#'
#' @noRd
#'
genm_is_minimum <- function(xyz, kmat) {
  rigid <- genm_rigid_basis(xyz)
  kmat_shifted <- kmat + mean(diag(kmat)) * tcrossprod(rigid)
  !is.null(tryCatch(chol(kmat_shifted), error = function(e) NULL))
}


#' Orthonormal basis of the rigid-body displacements
#'
#' Three translations and three infinitesimal rotations about the centroid. A
#' rotation about axis \eqn{u} moves the node at \eqn{r} by \eqn{u \times r}.
#'
#' @param xyz coordinates, a vector of length `3 * nsites`
#' @returns a `3 nsites x 6` matrix with orthonormal columns
#'
#' @noRd
#'
genm_rigid_basis <- function(xyz) {
  position <- matrix(xyz, nrow = 3)            # column k is node k
  r <- position - rowMeans(position)           # relative to the centroid
  nsites <- ncol(r)

  # each displacement as a 3 x nsites matrix, one column per node
  displacements <- list(
    move_x = matrix(c(1, 0, 0), nrow = 3, ncol = nsites),
    move_y = matrix(c(0, 1, 0), nrow = 3, ncol = nsites),
    move_z = matrix(c(0, 0, 1), nrow = 3, ncol = nsites),
    turn_x = rbind(0, -r[3, ], r[2, ]),          # (1, 0, 0) x r
    turn_y = rbind(r[3, ], 0, -r[1, ]),          # (0, 1, 0) x r
    turn_z = rbind(-r[2, ], r[1, ], 0)           # (0, 0, 1) x r
  )
  basis <- sapply(displacements, as.vector)     # 3N x 6
  qr.Q(qr(basis))                               # orthonormalised
}


#' Superpose coordinates onto a target (Kabsch)
#'
#' Applies to `xyz` the rotation and translation that minimise its RMSD to
#' `target`. Internal distances are unchanged.
#'
#' @param xyz,target coordinate vectors of length `3 * nsites`
#' @returns `xyz` superposed onto `target`
#'
#' @noRd
#'
genm_superpose <- function(xyz, target) {
  mobile <- matrix(xyz, nrow = 3)
  fixed <- matrix(target, nrow = 3)
  mobile_centroid <- rowMeans(mobile)
  fixed_centroid <- rowMeans(fixed)
  mobile <- mobile - mobile_centroid
  fixed <- fixed - fixed_centroid

  # Kabsch: the best rotation comes from the SVD of the 3 x 3 correlation
  # matrix; the sign correction keeps it a rotation rather than a reflection
  svd_corr <- svd(fixed %*% t(mobile))
  handedness <- sign(det(svd_corr$u %*% t(svd_corr$v)))
  rotation <- svd_corr$u %*% diag(c(1, 1, handedness)) %*% t(svd_corr$v)

  as.vector(rotation %*% mobile + fixed_centroid)
}


# Normal modes -----------------------------------------------------------------

#' Normal-mode analysis of a genm Hessian
#'
#' Same output and conventions as `calculate_enm_nma()` (eigenvalues ascending,
#' `mode = 1..nmodes`, [canonical_sign()]), but checks rather than assumes the
#' null space: exactly six eigenvalues with `|lambda| <= null_tol * max(lambda)`
#' and none below `-null_tol * max(lambda)`, else an error. A seventh null
#' direction is a network that is not rigid; a negative eigenvalue is not a
#' minimum.
#'
#' @param kmat the Hessian
#' @param null_tol relative threshold for a null eigenvalue
#'
#' @returns `list(mode, evalue, cmat, umat)`
#'
#' @noRd
#'
genm_nma <- function(kmat, null_tol = 1e-8) {
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
  list(
    mode = seq_along(evalue),
    evalue = evalue,
    cmat = umat %*% ((1 / evalue) * t(umat)), # pseudo-inverse of kmat
    umat = umat
  )
}


# Force constants -------------------------------------------------------------

#' The kij function of a genm
#'
#' Resolves `kij_<model>` with `match.fun()`, as `calculate_enm_graph()` does, and
#' checks that every name in `param$kij_par` is an argument it accepts.
#'
#' @param param the `param` list of a `"genm"` object
#' @returns the function
#'
#' @noRd
#'
genm_kij_fun <- function(param) {
  kij_fun <- tryCatch(match.fun(paste0("kij_", param$model)),
                      error = function(e) stop("unknown model: '", param$model, "'", call. = FALSE))
  extra <- param$kij_par
  if (length(extra) > 0) {
    if (is.null(names(extra)) || any(names(extra) == "")) {
      stop("parameters in ... must be named")
    }
    # what the user may set: the function's own parameters, not the ones penm passes
    accepted <- setdiff(names(formals(kij_fun)), c("dij", "sdij", "d_max", "..."))
    unknown <- setdiff(names(extra), accepted)
    if (length(unknown) > 0) {
      stop("kij_", param$model, " does not take parameter(s): ", paste(unknown, collapse = ", "),
           if (length(accepted) > 0) paste0(" (it takes: ", paste(accepted, collapse = ", "), ")"))
    }
  }
  kij_fun
}


#' Evaluate kij = k(lij)
#'
#' @param param the `param` list of a `"genm"` object
#' @param lij equilibrium lengths
#' @param sdij sequence separations
#'
#' @returns `kij`, checked to be finite and non-negative
#'
#' @noRd
#'
genm_kij <- function(param, lij, sdij) {
  kij_fun <- genm_kij_fun(param)
  # the kij_* functions name their first argument dij; here it receives lij
  kij <- do.call(kij_fun, c(list(lij, sdij = sdij, d_max = param$d_max), param$kij_par))
  if (length(kij) != length(lij)) {
    stop("kij_", param$model, " returned ", length(kij), " values for ", length(lij), " springs")
  }
  if (any(!is.finite(kij))) stop("kij_", param$model, " returned non-finite values")
  if (any(kij < 0)) {
    stop("kij_", param$model, " is negative for ", sum(kij < 0),
         " spring(s), e.g. at lij = ", signif(lij[kij < 0][1], 4), ": k must be >= 0")
  }
  kij
}


#' Warn when d_max_pairs truncates k
#'
#' Pairs of nodes farther apart than `d_max_pairs` in the pdb get no spring, so
#' they can never become contacts. That is harmless only if `k` is negligible
#' there. Compares `k(d_max_pairs)`, for a pair far apart in sequence, with the
#' largest `k` among springs more than three apart in sequence (beyond the
#' special-cased bonded neighbours of some models), and warns when the ratio
#' exceeds `ratio_max`. For power-law and gaussian models (pfanm, hnm0) there is
#' no distance where `k` is negligible, and the truncation then acts as a hard
#' cutoff.
#'
#' @param enm a `"genm"` object
#' @param ratio_max largest acceptable `k_far / k_near`
#'
#' @noRd
#'
genm_check_d_max_pairs <- function(enm, ratio_max = 0.01) {
  param <- enm$param
  far_in_sequence <- 1e6
  k_far <- genm_kij(param, param$d_max_pairs, sdij = far_in_sequence)
  k_near <- max(enm$springs$kij[enm$springs$sdij > 3])
  if (k_far > ratio_max * k_near) {
    warning("k(d_max_pairs = ", param$d_max_pairs, ") is ", signif(k_far / k_near, 2),
            " of the largest non-bonded k for model '", param$model, "'. ",
            "Pairs beyond d_max_pairs get no spring, so the truncation acts as a hard ",
            "cutoff and no contact can form beyond it.", call. = FALSE)
  }
  invisible(enm)
}


#' Check that a value is a single integer in a range
#'
#' @param x the value
#' @param name its name, for the error message
#' @param min,max the allowed range, inclusive
#'
#' @returns `x`, invisibly, or an error
#'
#' @noRd
#'
genm_check_integer <- function(x, name, min, max = Inf) {
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) &&
    x == round(x) && x >= min && x <= max
  if (!ok) {
    range <- if (is.finite(max)) paste0("in ", min, "..", max) else paste0(">= ", min)
    stop(name, " must be a single integer ", range)
  }
  invisible(x)
}
