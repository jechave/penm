# Generalized ENM (R/genm.R).
#
# Nothing here is frozen: every expectation is an invariant or a closed form.
# Where two quantities are compared, they are reached by different routes, so
# that the comparison is not equal by construction.

genm_mutate           <- penm:::genm_mutate
genm_allele_delta_lij <- penm:::genm_allele_delta_lij
genm_kij              <- penm:::genm_kij
genm_minimize         <- penm:::genm_minimize
genm_superpose_prot   <- penm:::genm_superpose_prot
genm_v_xyz            <- penm:::genm_v_xyz
genm_gradient         <- penm:::genm_gradient
genm_kmat             <- penm:::genm_kmat
set_enm_nma           <- penm:::set_enm_nma
dij_edge              <- penm:::dij_edge

load(test_path("fixtures", "pdb_2acy_A.rda"))

wt <- set_enm(pdb_2acy_A, node = "calpha", model = "ming_wall", d_max = 10.5,
              mut_model = "genm", d_max_graph = 14)
xyz_pdb <- penm:::calculate_enm_nodes(pdb_2acy_A, "ca")$xyz
nsites <- get_nsites(wt)

# length of every edge, at the protein's structure
edge_d <- function(prot) {
  dij_edge(get_xyz(prot), prot$graph$i, prot$graph$j)
}

# what the sequence determines
parameters <- function(prot) {
  list(sequence = prot$nodes$sequence, lij = prot$graph$lij, kij = prot$graph$kij)
}

# the edges of a site: their row in the graph, and the site at the other end
site_edges <- function(prot, site) {
  graph <- prot$graph
  row <- which(graph$i == site | graph$j == site)
  partner <- ifelse(graph$i[row] == site, graph$j[row], graph$i[row])
  tibble::tibble(row = row, partner = partner, sdij = graph$sdij[row],
                 l0ij = graph$l0ij[row], lij = graph$lij[row], kij = graph$kij[row])
}

# apply a path of mutations, given as a table with columns site and allele
mutate_path <- function(prot, path) {
  for (k in seq_len(nrow(path))) {
    prot <- genm_mutate(prot, site = path$site[k], allele = path$allele[k])
  }
  prot
}

superpose <- function(xyz, target) {
  all_coordinates <- seq_along(target)
  as.vector(bio3d::fit.xyz(fixed = target, mobile = xyz,
                           fixed.inds = all_coordinates, mobile.inds = all_coordinates))
}

# Two sites joined by an edge, so that mutating both changes a shared edge.
# site_b is the first site in contact with site_a (k > 0) that is not bonded to it.
site_a <- 80L
edges_a <- site_edges(wt, site_a)
site_b <- edges_a$partner[edges_a$kij > 0 & edges_a$sdij >= 2][1]

mut_a <- set_enm_nma(genm_mutate(wt, site_a, 3L))


# wild type -----------------------------------------------------------------

test_that("wild type: the pdb structure is the minimum, and V_min = 0", {
  expect_s3_class(wt, "prot")
  expect_identical(get_xyz(wt), xyz_pdb)
  expect_identical(enm_v_min(wt), 0)
  expect_identical(max(abs(genm_gradient(xyz_pdb, wt$graph, nsites))), 0)
  # the minimiser agrees: it stops at once
  p <- genm_minimize(wt)
  expect_equal(p$internal$minimization$iter, 0)
  expect_identical(get_xyz(p), xyz_pdb)
})

test_that("wild type: started off the structure, the minimiser returns to it", {
  # the test above exits at iteration 0 and so never exercises the minimiser
  seed <- xyz_pdb + 0.5 * sin(seq_along(xyz_pdb) * 1.7)
  expect_gt(genm_v_xyz(seed, wt$graph), 1)

  p <- wt
  p$nodes$xyz <- seed
  p <- genm_minimize(p)
  expect_gt(p$internal$minimization$iter, 0)
  expect_lt(enm_v_min(p), 1e-20)
  expect_lt(max(abs(edge_d(p) - wt$graph$lij)), 1e-9)
  expect_lt(max(abs(superpose(get_xyz(p), xyz_pdb) - xyz_pdb)), 1e-9)
})


# derivatives ---------------------------------------------------------------

test_that("the gradient matches finite differences of the energy", {
  # mutant parameters at the wild-type structure: far from stationary.
  # All edges, k = 0 included: the derivatives must hold for any set.
  edges <- mut_a$graph
  x <- xyz_pdb
  h <- 1e-5
  g <- genm_gradient(x, edges, nsites)
  g_fd <- vapply(seq_along(x), function(k) {
    e <- replace(numeric(length(x)), k, h)
    (genm_v_xyz(x + e, edges) - genm_v_xyz(x - e, edges)) / (2 * h)
  }, numeric(1))
  expect_gt(max(abs(g)), 0.1)
  expect_lt(max(abs(g - g_fd)), 1e-6)
})

test_that("the kmat, transverse term included, matches finite differences of the gradient", {
  # The only test of the transverse term: the six null modes and everything
  # else would survive a wrong sign on gij.
  edges <- mut_a$graph
  h <- 1e-5
  hess_fd <- function(x) {
    vapply(seq_along(x), function(k) {
      e <- replace(numeric(length(x)), k, h)
      (genm_gradient(x + e, edges, nsites) - genm_gradient(x - e, edges, nsites)) / (2 * h)
    }, numeric(length(x)))
  }

  # at the frustrated minimum, where the transverse term must be large enough
  # to be seen above the tolerance
  x <- get_xyz(mut_a)
  dij <- dij_edge(x, edges$i, edges$j)
  expect_gt(max(abs(edges$kij * (dij - edges$lij) / dij)), 0.1)
  expect_lt(max(abs(genm_kmat(x, edges, nsites) - hess_fd(x))), 1e-6)

  # and away from it: the formula holds at any conformation
  x <- xyz_pdb
  expect_lt(max(abs(genm_kmat(x, edges, nsites) - hess_fd(x))), 1e-6)
})


# rigidity ------------------------------------------------------------------

test_that("the kmat has exactly six null modes, for the wild type and for a mutant", {
  for (p in list(wt, mut_a)) {
    ev <- eigen(get_kmat(p), symmetric = TRUE, only.values = TRUE)$values
    thr <- 1e-8 * max(ev)
    expect_equal(sum(abs(ev) <= thr), 6)
    expect_gt(min(ev), -thr)
    expect_gt(sort(ev)[7], 1e3 * thr)   # a clear gap, not a marginal count
    expect_equal(get_nmodes(p), 3 * nsites - 6)
  }
})


# the sequence determines the parameters ------------------------------------

test_that("a protein's parameters are those its sequence implies, whatever the path", {
  # the same sequence reached by two paths: one mutates site_a twice and passes
  # through intermediate alleles, the other mutates each site once
  long_way <- mutate_path(wt, tibble::tibble(
    site   = c(site_a, site_b, site_a, 10, site_b),
    allele = c(3,      7,      5,      1,  2)
  ))
  short_way <- mutate_path(wt, tibble::tibble(
    site   = c(10, site_b, site_a),
    allele = c(1,  2,      5)
  ))
  expect_equal(long_way$nodes$sequence[c(site_a, site_b, 10)], c(5L, 2L, 1L))
  expect_equal(sum(long_way$nodes$sequence != 0), 3)
  expect_identical(parameters(long_way), parameters(short_way))
})

test_that("state function: two orders, through different intermediates, reach the same protein", {
  expect_true(any(wt$graph$i == min(site_a, site_b) & wt$graph$j == max(site_a, site_b) &
                    wt$graph$kij > 0))

  # a then b, and b then a, each step starting from the previous structure
  p_ab <- genm_mutate(mut_a, site_b, 7L)
  p_b <- genm_mutate(wt, site_b, 7L)
  p_ba <- genm_mutate(p_b, site_a, 3L)

  # the routes genuinely differ
  expect_gt(abs(enm_v_min(mut_a) - enm_v_min(p_b)), 0.1)
  expect_gt(max(abs(get_xyz(mut_a) - get_xyz(p_b))), 0.01)
  expect_gt(p_ab$internal$minimization$iter, 0)
  expect_gt(p_ba$internal$minimization$iter, 0)

  # same parameters, bit for bit; same minimum, to minimiser precision
  expect_identical(parameters(p_ab), parameters(p_ba))
  expect_lt(abs(enm_v_min(p_ab) - enm_v_min(p_ba)) / enm_v_min(p_ab), 1e-12)
  expect_lt(max(abs(edge_d(p_ab) - edge_d(p_ba))), 1e-9)
})

test_that("reversibility: mutating back restores the protein, from the wild type and from a frustrated reference", {
  # from the wild type
  back <- genm_mutate(mut_a, site_a, 0L)
  expect_identical(parameters(back), parameters(wt))
  expect_gt(enm_v_min(mut_a), 0.1)
  expect_gt(back$internal$minimization$iter, 0)
  expect_lt(enm_v_min(back), 1e-20)
  expect_lt(max(abs(edge_d(back) - edge_d(wt))), 1e-9)

  # from a frustrated reference: the mutant at site_b
  ref <- genm_mutate(wt, site_b, 7L)
  fwd <- genm_mutate(ref, site_a, 3L)
  back <- genm_mutate(fwd, site_a, 0L)
  expect_identical(parameters(back), parameters(ref))
  expect_gt(enm_v_min(ref), 0.1)
  expect_gt(abs(enm_v_min(fwd) - enm_v_min(ref)), 0.1)
  expect_lt(abs(enm_v_min(back) - enm_v_min(ref)) / enm_v_min(ref), 1e-12)
  expect_lt(max(abs(edge_d(back) - edge_d(ref))), 1e-9)
})

test_that("a path out and back, from an arbitrary sequence, returns to it exactly", {
  # Start away from the wild type, with nonzero alleles on both ends of some
  # edges: returning to all zeros would be the degenerate case, where every
  # delta vanishes and errors in combining them cannot show.
  start <- mutate_path(wt, tibble::tibble(
    site   = c(site_b, 10, 40),
    allele = c(7,      1,  9)
  ))

  # out: re-mutate a site that is already mutated, mutate site_a (which shares
  # an edge with site_b) twice, and add new sites
  p <- mutate_path(start, tibble::tibble(
    site   = c(site_a, site_b, 5, site_a, 40, 70),
    allele = c(3,      2,      4, 6,      1,  8)
  ))
  expect_gt(max(abs(p$graph$lij - start$graph$lij)), 0.1)

  # back, in a different order from the way out
  for (s in rev(which(p$nodes$sequence != start$nodes$sequence))) {
    p <- genm_mutate(p, site = s, allele = start$nodes$sequence[s])
  }
  expect_identical(parameters(p), parameters(start))
})


# the mutational process ----------------------------------------------------

test_that("allele 0 changes nothing, and other alleles change exactly the edges they should", {
  expect_identical(genm_allele_delta_lij(wt, site_a, 0L), numeric(nsites))

  # allele 3 at site_a changes the edges of site_a with sdij >= 2, and no other
  changed <- which(mut_a$graph$lij != wt$graph$lij)
  expect_setequal(changed, edges_a$row[edges_a$sdij >= 2])
  expect_gt(length(changed), 10)

  # mut_sd_min = 1 perturbs the i,i+1 bonds as well
  wt_1 <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
                  mut_model = "genm", d_max_graph = 14, mut_sd_min = 1)
  changed_1 <- which(genm_mutate(wt_1, site_a, 3L)$graph$lij != wt_1$graph$lij)
  expect_setequal(changed_1, site_edges(wt_1, site_a)$row)
})

test_that("an allele's draw depends on (ensemble, site, allele), and leaves the caller's RNG alone", {
  d1 <- genm_allele_delta_lij(wt, site_a, 3L)
  expect_identical(genm_allele_delta_lij(wt, site_a, 3L), d1)
  expect_false(isTRUE(all.equal(genm_allele_delta_lij(wt, site_a, 4L), d1)))
  expect_false(isTRUE(all.equal(genm_allele_delta_lij(wt, site_b, 3L), d1)))
  wt_2 <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
                  mut_model = "genm", d_max_graph = 14, ensemble = 2L)
  expect_false(isTRUE(all.equal(genm_allele_delta_lij(wt_2, site_a, 3L), d1)))
  # sd of the draws is mut_dl_sigma (0.3), loosely: nsites draws
  expect_gt(sd(d1), 0.2)
  expect_lt(sd(d1), 0.4)

  set.seed(42); x <- runif(3)
  set.seed(42); invisible(genm_mutate(wt, site_a, 3L)); y <- runif(3)
  expect_identical(x, y)
})

test_that("an allele's change to an edge does not depend on which other edges exist", {
  # the change allele 3 at site_a makes to each of its edges, labelled by the
  # site at the edge's other end
  change_by_partner <- function(p) {
    edges <- site_edges(genm_mutate(p, site_a, 3L), site_a)
    tibble::tibble(partner = edges$partner, sdij = edges$sdij,
                   change = edges$lij - edges$l0ij)
  }
  narrow <- change_by_partner(wt)
  perturbed <- narrow$sdij >= 2
  # the shared edges are perturbed, so equal changes are not 0 == 0
  expect_gt(sum(narrow$change[perturbed] != 0), 10)

  # more edges: site_a gains partners, and the ones it had keep their change
  wide <- change_by_partner(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
                                    mut_model = "genm", d_max_graph = 16))
  expect_gt(nrow(wide), nrow(narrow))
  same_edge_in_wide <- match(narrow$partner, wide$partner)
  expect_false(anyNA(same_edge_in_wide))
  expect_identical(wide$change[same_edge_in_wide], narrow$change)

  # mut_sd_min = 1 perturbs the i,i+1 bonds too; the edges both perturb keep
  # their change
  all_perturbed <- change_by_partner(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
                                             mut_model = "genm", d_max_graph = 14, mut_sd_min = 1))
  expect_identical(all_perturbed$partner, narrow$partner)
  expect_identical(all_perturbed$change[perturbed], narrow$change[perturbed])
})

test_that("genm_mutate validates its input", {
  expect_error(genm_mutate(mut_a, site_a, 3L), "already has allele")
  expect_error(genm_mutate(wt, site_a, 0L), "already has allele")
  expect_error(genm_mutate(wt, site_a, 10L), "allele must be")
  expect_error(genm_mutate(wt, site_a, -1L), "allele must be")
  expect_error(genm_mutate(wt, site_a, 1.5), "allele must be")
  expect_error(genm_mutate(wt, nsites + 1, 1L), "site must be")
  expect_error(genm_mutate(wt, 1.5, 1L), "site must be")
  # a draw so wide that some length would become negative
  wide <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
                  mut_model = "genm", d_max_graph = 14, mut_dl_sigma = 50)
  expect_error(genm_mutate(wide, site_a, 3L), "<= 0")
  # a prot built for lfenm
  lfenm_wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5)
  expect_error(genm_mutate(lfenm_wt, site_a, 3L), "mut_model")
  # and get_mutant_site does not mutate a genm prot yet
  expect_error(get_mutant_site(wt, site_a, 3L), "does not support")
})


# parameters ----------------------------------------------------------------

test_that("kij follows lij: a mutation that moves an edge across the cutoff changes its k", {
  # search the alleles of site_a for one whose draw carries an edge across
  # d_max = 10.5 in either direction
  crossed <- NULL
  for (a in 1:9) {
    m <- genm_mutate(wt, site_a, a)
    flip <- which((m$graph$lij[edges_a$row] <= 10.5) != (edges_a$lij <= 10.5) & edges_a$sdij >= 2)
    if (length(flip) > 0) { crossed <- list(m = m, rows = edges_a$row[flip]); break }
  }
  expect_false(is.null(crossed))
  m <- crossed$m
  r <- crossed$rows
  expect_true(all(m$graph$kij[r] != wt$graph$kij[r]))
  expect_equal(m$graph$kij[r], ifelse(m$graph$lij[r] <= 10.5, 4.5, 0))
  # and in general, kij is k(lij) everywhere
  expect_equal(m$graph$kij, penm:::kij_ming_wall(m$graph$lij, m$graph$sdij, d_max = 10.5))
})

test_that("a negative k is rejected", {
  # kij_hnm is negative below 2390/860 = 2.78 A
  param <- list(model = "hnm", d_max = 10.5, kij_par = list())
  expect_error(genm_kij(param, c(5, 2.5), c(5, 5)), "k must be >= 0")
  expect_no_error(genm_kij(param, c(5, 3), c(5, 5)))
})

test_that("set_enm validates the model, its parameters and the mutational process, for genm", {
  b <- function(...) set_enm(pdb_2acy_A, node = "ca", d_max = 10.5, mut_model = "genm", ...)
  expect_error(b(model = "no_such_model", d_max_graph = 14), "kij_no_such_model")
  expect_error(b(model = "ming_wall", d_max_graph = 14, wdth = 1), "does not take parameter")
  expect_error(b(model = "ming_wall", d_max_graph = 14, 1), "must be named")
  expect_error(b(model = "ming_wall_smooth", d_max_graph = 14), "\"d_max_width\" is missing")
  expect_error(b(model = "ming_wall", d_max_graph = 9), "must not be smaller")
  expect_identical(b(model = "ming_wall_smooth", d_max_graph = 14, d_max_width = 1)$param$kij_par,
                   list(d_max_width = 1))
  expect_error(b(model = "ming_wall", d_max_graph = 14, ensemble = NA), "ensemble")
  expect_error(b(model = "ming_wall", d_max_graph = 14, n_alleles = 1), "n_alleles")
  expect_error(b(model = "ming_wall", d_max_graph = 14, mut_dl_sigma = 0), "mut_dl_sigma")
  expect_error(b(model = "ming_wall", d_max_graph = 14, mut_sd_min = 0), "mut_sd_min")
})

test_that("there is an edge for every pair within d_max_graph and every i,i+1 pair", {
  x <- matrix(xyz_pdb, 3)
  d <- as.matrix(dist(t(x)))
  bonded <- abs(outer(get_pdb_site(wt), get_pdb_site(wt), "-")) == 1
  expect_equal(nrow(wt$graph), sum((d <= 14 | bonded)[upper.tri(d)]))
  expect_true(all(wt$graph$i < wt$graph$j))
  expect_false(is.unsorted(wt$graph$i * (nsites + 1) + wt$graph$j))
  # l0ij is the pdb distance, and at allele 0 everywhere lij = l0ij
  expect_equal(wt$graph$l0ij, d[cbind(wt$graph$i, wt$graph$j)])
  expect_identical(wt$graph$lij, wt$graph$l0ij)
  expect_identical(wt$nodes$sequence, integer(nsites))
  # every i,i+1 pair is present, with sdij = 1
  ij <- which(bonded & upper.tri(bonded), arr.ind = TRUE)
  key <- paste(wt$graph$i, wt$graph$j)
  expect_true(all(paste(ij[, 1], ij[, 2]) %in% key))
  expect_true(all(wt$graph$sdij[key %in% paste(ij[, 1], ij[, 2])] == 1))
})

test_that("d_max_graph warns when the truncation is a hard cutoff", {
  b <- function(...) set_enm(pdb_2acy_A, node = "ca", d_max = 10.5, mut_model = "genm", d_max_graph = 14, ...)
  expect_warning(b(model = "pfanm"), "hard cutoff")
  expect_warning(b(model = "hnm0"), "hard cutoff")
  expect_no_warning(b(model = "ming_wall"))
  expect_no_warning(b(model = "ming_wall_smooth", d_max_width = 1))
})


# the protein ---------------------------------------------------------------

test_that("the minimum is superposed onto the seed", {
  # the mutant moved from its seed, the wild type's structure ...
  expect_gt(max(abs(get_xyz(mut_a) - get_xyz(wt))), 0.01)
  # ... and is already superposed onto it: superposing again is a no-op
  expect_lt(max(abs(superpose(get_xyz(mut_a), get_xyz(wt)) - get_xyz(mut_a))), 1e-10)
})

test_that("dij in the graph are the edge lengths of the protein's structure", {
  expect_identical(mut_a$graph$dij, edge_d(mut_a))
  expect_gt(max(abs(mut_a$graph$dij - wt$graph$dij)), 0.01)
})

test_that("genm_superpose_prot rotates the whole protein, kmat and modes included", {
  # a known rotation and translation of the mutant's own structure
  th <- 0.7
  rot <- matrix(c(cos(th), sin(th), 0, -sin(th), cos(th), 0, 0, 0, 1), 3)
  moved <- as.vector(rot %*% matrix(get_xyz(mut_a), 3) + c(5, -3, 2))

  q <- genm_superpose_prot(mut_a, target = moved)
  expect_lt(max(abs(get_xyz(q) - moved)), 1e-10)

  # The kmat and the covariance must be the original ones, rotated: with the
  # coordinates ordered x1 y1 z1 x2 ..., the rotation of the whole protein is
  # block-diagonal, one copy of rot per node.
  rot_all <- kronecker(diag(nsites), rot)
  k <- get_kmat(mut_a)
  expect_lt(max(abs(get_kmat(q) - rot_all %*% k %*% t(rot_all))), 1e-9 * max(abs(k)))
  cmat <- get_cmat(mut_a)
  expect_lt(max(abs(get_cmat(q) - rot_all %*% cmat %*% t(rot_all))), 1e-9 * max(abs(cmat)))
  expect_equal(get_evalue(q), get_evalue(mut_a))

  # what does not depend on orientation is unchanged
  expect_identical(parameters(q), parameters(mut_a))
  expect_identical(q$graph$dij, edge_d(q))
  expect_lt(max(abs(q$graph$dij - mut_a$graph$dij)), 1e-10)
  expect_lt(abs(enm_v_min(q) - enm_v_min(mut_a)) / enm_v_min(mut_a), 1e-12)

  # a protein without modes stays without modes
  no_modes <- genm_mutate(wt, site_a, 3L)
  expect_identical(genm_superpose_prot(no_modes, target = moved)$nma, NA)

  expect_error(genm_superpose_prot(mut_a, target = moved[-1]), "length")
})

test_that("the existing analysis getters work on a genm prot", {
  expect_equal(length(get_msf_site(mut_a)), nsites)
  expect_true(all(is.finite(get_msf_site(mut_a)) & get_msf_site(mut_a) > 0))
  # cmat is the pseudo-inverse of kmat
  k <- get_kmat(mut_a)
  expect_lt(max(abs(k %*% get_cmat(mut_a) %*% k - k)), 1e-8 * max(abs(k)))
  # the mutation changes the spectrum, unlike lfenm
  expect_gt(abs(ddg_tds(wt, mut_a)), 0)
  # eij belongs to lfenm, and a genm prot has none
  expect_error(penm:::get_eij(wt), "no eij")
})

test_that("genm_mutate leaves no modes; set_enm_nma adds them and changes nothing else", {
  p <- genm_mutate(wt, site_a, 3L)
  expect_identical(p$nma, NA)
  expect_error(get_evalue(p))
  expect_error(get_cmat(p))

  q <- set_enm_nma(p)
  expect_equal(get_nmodes(q), 3 * nsites - 6)
  expect_identical(q[names(q) != "nma"], p[names(p) != "nma"])
  expect_identical(class(q), class(p))
  # recomputing on a protein that has modes gives the same modes
  expect_identical(set_enm_nma(q), q)

  # a mutant of a protein with modes does not keep them: they would be stale
  expect_identical(genm_mutate(mut_a, site_b, 7L)$nma, NA)
})
