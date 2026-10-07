# Generalized ENM (R/genm.R).
#
# Nothing here is frozen: every expectation is an invariant or a closed form, so
# there are no fixtures to regenerate. Where two quantities are compared, the
# test also asserts that they were reached by different routes, so that the
# comparison is not equal by construction.

build_enm_from_pdb    <- penm:::build_enm_from_pdb
genm_mutate           <- penm:::genm_mutate
genm_lij              <- penm:::genm_lij
genm_site_dl          <- penm:::genm_site_dl
prot_from_enm         <- penm:::prot_from_enm
genm_site_springs     <- penm:::genm_site_springs
genm_add_nma          <- penm:::genm_add_nma
genm_v_min            <- penm:::genm_v_min
genm_energy           <- penm:::genm_energy
genm_gradient         <- penm:::genm_gradient
genm_hessian          <- penm:::genm_hessian
genm_geometry         <- penm:::genm_geometry
genm_superpose        <- penm:::genm_superpose
genm_superpose_prot   <- penm:::genm_superpose_prot
genm_kij              <- penm:::genm_kij

load(test_path("fixtures", "pdb_2acy_A.rda"))

enm <- build_enm_from_pdb(pdb_2acy_A, node = "calpha", model = "ming_wall",
                          d_max = 10.5, d_max_pairs = 14)
xyz_pdb <- penm:::calculate_enm_nodes(pdb_2acy_A, "ca")$xyz
wt <- genm_add_nma(prot_from_enm(enm, xyz_pdb))
nsites <- get_nsites(wt)

# length of every spring, at the protein's minimum
spring_d <- function(prot) {
  genm_geometry(get_xyz(prot), prot$enm$springs$i, prot$enm$springs$j)$dij
}

# apply a path of mutations, given as a table with columns site and allele
mutate_path <- function(enm, path) {
  for (k in seq_len(nrow(path))) {
    enm <- genm_mutate(enm, site = path$site[k], allele = path$allele[k])
  }
  enm
}

# Two sites joined by a spring, so that mutating both changes a shared spring.
# site_b is the first site in contact with site_a (k > 0) that is not bonded to it.
site_a <- 80L
springs_a <- genm_site_springs(enm, site_a)
contacts_a <- springs_a[springs_a$kij > 0 & springs_a$sdij >= 2, ]
site_b <- if (contacts_a$i[1] == site_a) contacts_a$j[1] else contacts_a$i[1]

mut_a_enm <- genm_mutate(enm, site_a, 3L)
mut_a <- genm_add_nma(prot_from_enm(mut_a_enm, get_xyz(wt)))


# wild type -----------------------------------------------------------------

test_that("wild type: the minimiser does not move the pdb structure, and V_min = 0", {
  expect_equal(wt$minimization$iter, 0)
  expect_identical(get_xyz(wt), xyz_pdb)
  expect_identical(genm_v_min(wt), 0)
})

test_that("wild type: started off the structure, the minimiser returns to it", {
  # the test above exits at iteration 0 and so never exercises the minimiser
  seed <- xyz_pdb + 0.5 * sin(seq_along(xyz_pdb) * 1.7)
  expect_gt(genm_energy(seed, enm$springs), 1)

  p <- prot_from_enm(enm, seed)
  expect_gt(p$minimization$iter, 0)
  expect_lt(genm_v_min(p), 1e-20)
  expect_lt(max(abs(spring_d(p) - enm$springs$lij)), 1e-9)
  expect_lt(max(abs(genm_superpose(get_xyz(p), xyz_pdb) - xyz_pdb)), 1e-9)
})


# derivatives ---------------------------------------------------------------

test_that("the gradient matches finite differences of the energy", {
  # mutant parameters at the wild-type structure: far from stationary.
  # All springs, k = 0 included: the derivatives must hold for any set.
  sp <- mut_a_enm$springs
  x <- xyz_pdb
  h <- 1e-5
  g <- genm_gradient(x, sp, nsites)
  g_fd <- vapply(seq_along(x), function(k) {
    e <- replace(numeric(length(x)), k, h)
    (genm_energy(x + e, sp) - genm_energy(x - e, sp)) / (2 * h)
  }, numeric(1))
  expect_gt(max(abs(g)), 0.1)
  expect_lt(max(abs(g - g_fd)), 1e-6)
})

test_that("the Hessian, transverse term included, matches finite differences of the gradient", {
  # The only test of the transverse term: the six null modes and everything
  # else would survive a wrong sign on gij.
  sp <- mut_a_enm$springs
  h <- 1e-5
  hess_fd <- function(x) {
    vapply(seq_along(x), function(k) {
      e <- replace(numeric(length(x)), k, h)
      (genm_gradient(x + e, sp, nsites) - genm_gradient(x - e, sp, nsites)) / (2 * h)
    }, numeric(length(x)))
  }

  # at the frustrated minimum, where the transverse term must be large enough
  # to be seen above the tolerance
  x <- get_xyz(mut_a)
  geo <- genm_geometry(x, sp$i, sp$j)
  expect_gt(max(abs(sp$kij * (geo$dij - sp$lij) / geo$dij)), 0.1)
  expect_lt(max(abs(genm_hessian(x, sp, nsites) - hess_fd(x))), 1e-6)

  # and away from it: the formula holds at any conformation
  x <- xyz_pdb
  expect_lt(max(abs(genm_hessian(x, sp, nsites) - hess_fd(x))), 1e-6)
})


# rigidity ------------------------------------------------------------------

test_that("the Hessian has exactly six null modes, for the wild type and for a mutant", {
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

test_that("an enm's parameters are those its sequence implies, whatever the path", {
  # a path that mutates a site twice and two sites that share a spring
  e <- mutate_path(enm, tibble::tibble(
    site   = c(site_a, site_b, site_a, 10, site_b),
    allele = c(3,      7,      5,      1,  2)
  ))
  expect_equal(e$sequence[c(site_a, site_b, 10)], c(5L, 2L, 1L))
  expect_equal(sum(e$sequence != 0), 3)
  # lij recomputed from l0ij and the sequence alone matches, bit for bit
  expect_identical(e$springs$lij, genm_lij(e, seq_len(nrow(e$springs))))
  expect_identical(e$springs$kij, genm_kij(e$param, e$springs$lij, e$springs$sdij))
})

test_that("state function: two orders, through different intermediates, reach the same protein", {
  expect_true(any(enm$springs$i == min(site_a, site_b) & enm$springs$j == max(site_a, site_b) &
                    enm$springs$kij > 0))

  # a then b, each step seeded from the previous structure
  e_ab <- genm_mutate(mut_a_enm, site_b, 7L)
  p_ab <- prot_from_enm(e_ab, get_xyz(mut_a))

  # b then a
  e_b <- genm_mutate(enm, site_b, 7L)
  p_b <- prot_from_enm(e_b, get_xyz(wt))
  e_ba <- genm_mutate(e_b, site_a, 3L)
  p_ba <- prot_from_enm(e_ba, get_xyz(p_b))

  # the routes genuinely differ
  expect_gt(abs(genm_v_min(mut_a) - genm_v_min(p_b)), 0.1)
  expect_gt(max(abs(get_xyz(mut_a) - get_xyz(p_b))), 0.01)
  expect_gt(p_ab$minimization$iter, 0)
  expect_gt(p_ba$minimization$iter, 0)

  # same parameters, bit for bit; same minimum, to minimiser precision
  expect_identical(e_ab, e_ba)
  expect_lt(abs(genm_v_min(p_ab) - genm_v_min(p_ba)) / genm_v_min(p_ab), 1e-12)
  expect_lt(max(abs(spring_d(p_ab) - spring_d(p_ba))), 1e-9)
})

test_that("reversibility: mutating back restores the protein, from the wild type and from a frustrated reference", {
  # from the wild type
  back_enm <- genm_mutate(mut_a_enm, site_a, 0L)
  back <- prot_from_enm(back_enm, get_xyz(mut_a))
  expect_identical(back_enm, enm)
  expect_gt(genm_v_min(mut_a), 0.1)
  expect_gt(back$minimization$iter, 0)
  expect_lt(genm_v_min(back), 1e-20)
  expect_lt(max(abs(spring_d(back) - spring_d(wt))), 1e-9)

  # from a frustrated reference: the mutant at site_b
  ref_enm <- genm_mutate(enm, site_b, 7L)
  ref <- prot_from_enm(ref_enm, get_xyz(wt))
  fwd_enm <- genm_mutate(ref_enm, site_a, 3L)
  fwd <- prot_from_enm(fwd_enm, get_xyz(ref))
  back_enm <- genm_mutate(fwd_enm, site_a, 0L)
  back <- prot_from_enm(back_enm, get_xyz(fwd))
  expect_identical(back_enm, ref_enm)
  expect_gt(genm_v_min(ref), 0.1)
  expect_gt(abs(genm_v_min(fwd) - genm_v_min(ref)), 0.1)
  expect_lt(abs(genm_v_min(back) - genm_v_min(ref)) / genm_v_min(ref), 1e-12)
  expect_lt(max(abs(spring_d(back) - spring_d(ref))), 1e-9)
})

test_that("a path out and back, from an arbitrary sequence, returns to it exactly", {
  # Start away from the wild type, with nonzero alleles on both ends of some
  # springs: returning to all zeros would be the degenerate case, where every
  # delta vanishes and errors in combining them cannot show.
  start <- mutate_path(enm, tibble::tibble(
    site   = c(site_b, 10, 40),
    allele = c(7,      1,  9)
  ))

  # out: re-mutate a site that is already mutated, mutate site_a (which shares a
  # spring with site_b) twice, and add new sites
  e <- mutate_path(start, tibble::tibble(
    site   = c(site_a, site_b, 5, site_a, 40, 70),
    allele = c(3,      2,      4, 6,      1,  8)
  ))
  expect_gt(max(abs(e$springs$lij - start$springs$lij)), 0.1)

  # back, in a different order from the way out
  for (s in rev(which(e$sequence != start$sequence))) {
    e <- genm_mutate(e, site = s, allele = start$sequence[s])
  }
  expect_identical(e, start)
})


# the mutational process ----------------------------------------------------

test_that("allele 0 changes nothing, and other alleles change exactly the springs they should", {
  expect_identical(genm_site_dl(enm, site_a, 0L), numeric(nrow(genm_site_springs(enm, site_a))))

  sp <- genm_site_springs(enm, site_a)
  dl <- genm_site_dl(enm, site_a, 3L)
  expect_true(all(dl[sp$sdij < 2] == 0))
  expect_true(all(dl[sp$sdij >= 2] != 0))
  # and only those: every other spring is untouched
  changed <- which(mut_a_enm$springs$lij != enm$springs$lij)
  expect_setequal(changed, sp$pair[sp$sdij >= 2])

  # mut_sd_min = 1 perturbs the i,i+1 bonds as well
  e1 <- build_enm_from_pdb(pdb_2acy_A, node = "ca", model = "ming_wall",
                           d_max = 10.5, d_max_pairs = 14, mut_sd_min = 1)
  expect_true(all(genm_site_dl(e1, site_a, 3L) != 0))
})

test_that("an allele's draw depends on (ensemble, site, allele), and leaves the caller's RNG alone", {
  d1 <- genm_site_dl(enm, site_a, 3L)
  expect_identical(genm_site_dl(enm, site_a, 3L), d1)
  expect_false(isTRUE(all.equal(genm_site_dl(enm, site_a, 4L), d1)))
  e2 <- build_enm_from_pdb(pdb_2acy_A, node = "ca", model = "ming_wall",
                           d_max = 10.5, d_max_pairs = 14, ensemble = 2L)
  expect_false(isTRUE(all.equal(genm_site_dl(e2, site_a, 3L), d1)))
  # sd of the draws is mut_dl_sigma (0.3), loosely: about 40 draws
  expect_gt(sd(d1[d1 != 0]), 0.15)
  expect_lt(sd(d1[d1 != 0]), 0.45)

  set.seed(42); x <- runif(3)
  set.seed(42); invisible(genm_mutate(enm, site_a, 3L)); y <- runif(3)
  expect_identical(x, y)
})

test_that("an allele's change to a spring does not depend on which other springs exist", {
  # the change allele 3 at site_a makes to each of its springs, labelled by the
  # site at the spring's other end
  dl_by_partner <- function(e) {
    springs <- genm_site_springs(e, site_a)
    partner <- ifelse(springs$i == site_a, springs$j, springs$i)
    dl <- genm_site_dl(e, site_a, 3L)
    tibble::tibble(partner = partner, sdij = springs$sdij, dl = dl)
  }
  narrow <- dl_by_partner(enm)
  perturbed <- narrow$sdij >= 2
  # the shared springs are perturbed, so equal changes are not 0 == 0
  expect_gt(sum(narrow$dl[perturbed] != 0), 10)

  # more springs: site_a gains partners, and the ones it had keep their change
  wide <- dl_by_partner(build_enm_from_pdb(pdb_2acy_A, node = "ca", model = "ming_wall",
                                           d_max = 10.5, d_max_pairs = 16))
  expect_gt(nrow(wide), nrow(narrow))
  same_spring_in_wide <- match(narrow$partner, wide$partner)
  expect_false(anyNA(same_spring_in_wide))
  expect_identical(wide$dl[same_spring_in_wide], narrow$dl)

  # mut_sd_min = 1 perturbs the i,i+1 bonds too; the springs both perturb keep
  # their change
  all_perturbed <- dl_by_partner(build_enm_from_pdb(pdb_2acy_A, node = "ca", model = "ming_wall",
                                                    d_max = 10.5, d_max_pairs = 14,
                                                    mut_sd_min = 1))
  expect_identical(all_perturbed$partner, narrow$partner)
  expect_identical(all_perturbed$dl[perturbed], narrow$dl[perturbed])
})

test_that("genm_mutate validates its input", {
  expect_error(genm_mutate(mut_a_enm, site_a, 3L), "already has allele")
  expect_error(genm_mutate(enm, site_a, 0L), "already has allele")
  expect_error(genm_mutate(enm, site_a, 10L), "allele must be")
  expect_error(genm_mutate(enm, site_a, -1L), "allele must be")
  expect_error(genm_mutate(enm, site_a, 1.5), "allele must be")
  expect_error(genm_mutate(enm, nsites + 1, 1L), "site must be")
  expect_error(genm_mutate(enm, 1.5, 1L), "site must be")
  # a draw so wide that some length would become negative
  wide <- build_enm_from_pdb(pdb_2acy_A, node = "ca", model = "ming_wall",
                             d_max = 10.5, d_max_pairs = 14, mut_dl_sigma = 50)
  expect_error(genm_mutate(wide, site_a, 3L), "<= 0")
})


# parameters ----------------------------------------------------------------

test_that("kij follows lij: a mutation that moves a spring across the cutoff changes its k", {
  # search the alleles of site_a for one whose draw carries a spring across
  # d_max = 10.5 in either direction
  sp <- genm_site_springs(enm, site_a)
  crossed <- NULL
  for (a in 1:9) {
    m <- genm_mutate(enm, site_a, a)
    flip <- which((m$springs$lij[sp$pair] <= 10.5) != (sp$lij <= 10.5) & sp$sdij >= 2)
    if (length(flip) > 0) { crossed <- list(m = m, rows = sp$pair[flip]); break }
  }
  expect_false(is.null(crossed))
  m <- crossed$m
  r <- crossed$rows
  expect_true(all(m$springs$kij[r] != enm$springs$kij[r]))
  expect_equal(m$springs$kij[r], ifelse(m$springs$lij[r] <= 10.5, 4.5, 0))
  # and in general, kij is k(lij) everywhere
  expect_equal(m$springs$kij, penm:::kij_ming_wall(m$springs$lij, m$springs$sdij, d_max = 10.5))
})

test_that("a negative k is rejected", {
  # kij_hnm is negative below 2390/860 = 2.78 A
  param <- list(model = "hnm", d_max = 10.5, kij_par = list())
  expect_error(genm_kij(param, c(5, 2.5), c(5, 5)), "k must be >= 0")
  expect_no_error(genm_kij(param, c(5, 3), c(5, 5)))
})

test_that("build_enm_from_pdb validates the model, its parameters and the mutational process", {
  b <- function(...) build_enm_from_pdb(pdb_2acy_A, node = "ca", d_max = 10.5, d_max_pairs = 14, ...)
  expect_error(b(model = "no_such_model"), "unknown model")
  expect_error(b(model = "ming_wall", wdth = 1), "does not take parameter")
  expect_error(b(model = "ming_wall", 1), "must be named")
  expect_error(b(model = "ming_wall_smooth"), "\"w\" is missing")
  expect_error(build_enm_from_pdb(pdb_2acy_A, node = "ca", model = "ming_wall",
                                  d_max = 10.5, d_max_pairs = 9), "must not be smaller")
  expect_identical(b(model = "ming_wall_smooth", w = 1)$param$kij_par, list(w = 1))
  expect_error(b(model = "ming_wall", ensemble = NA), "ensemble")
  expect_error(b(model = "ming_wall", n_alleles = 1), "n_alleles")
  expect_error(b(model = "ming_wall", mut_dl_sigma = 0), "mut_dl_sigma")
  expect_error(b(model = "ming_wall", mut_sd_min = 0), "mut_sd_min")
})

test_that("there is a spring for every pair within d_max_pairs and every i,i+1 pair", {
  x <- matrix(xyz_pdb, 3)
  d <- as.matrix(dist(t(x)))
  bonded <- abs(outer(get_pdb_site(wt), get_pdb_site(wt), "-")) == 1
  expect_equal(nrow(enm$springs), sum((d <= 14 | bonded)[upper.tri(d)]))
  expect_true(all(enm$springs$i < enm$springs$j))
  expect_false(is.unsorted(enm$springs$i * (nsites + 1) + enm$springs$j))
  # l0ij is the pdb distance, and at allele 0 everywhere lij = l0ij
  expect_equal(enm$springs$l0ij, d[cbind(enm$springs$i, enm$springs$j)])
  expect_identical(enm$springs$lij, enm$springs$l0ij)
  expect_identical(enm$sequence, integer(nsites))
  # every i,i+1 pair is present, with sdij = 1
  ij <- which(bonded & upper.tri(bonded), arr.ind = TRUE)
  key <- paste(enm$springs$i, enm$springs$j)
  expect_true(all(paste(ij[, 1], ij[, 2]) %in% key))
  expect_true(all(enm$springs$sdij[key %in% paste(ij[, 1], ij[, 2])] == 1))
})

test_that("d_max_pairs warns when the truncation is a hard cutoff", {
  b <- function(...) build_enm_from_pdb(pdb_2acy_A, node = "ca", d_max = 10.5, d_max_pairs = 14, ...)
  expect_warning(b(model = "pfanm"), "hard cutoff")
  expect_warning(b(model = "hnm0"), "hard cutoff")
  expect_no_warning(b(model = "ming_wall"))
  expect_no_warning(b(model = "ming_wall_smooth", w = 1))
})


# the protein ---------------------------------------------------------------

test_that("the minimum is optimally superposed onto the seed", {
  # closed form: a known rotation and translation are undone exactly
  th <- 0.7
  rot <- matrix(c(cos(th), sin(th), 0, -sin(th), cos(th), 0, 0, 0, 1), 3)
  moved <- as.vector(rot %*% matrix(xyz_pdb, 3) + c(5, -3, 2))
  expect_lt(max(abs(genm_superpose(moved, xyz_pdb) - xyz_pdb)), 1e-10)

  # the mutant is already superposed onto its seed: superposing again is a no-op
  expect_lt(max(abs(genm_superpose(get_xyz(mut_a), get_xyz(wt)) - get_xyz(mut_a))), 1e-10)
})

test_that("genm_superpose_prot rotates the whole protein, Hessian and modes included", {
  # a known rotation and translation of the mutant's own structure
  th <- 0.7
  rot <- matrix(c(cos(th), sin(th), 0, -sin(th), cos(th), 0, 0, 0, 1), 3)
  moved <- as.vector(rot %*% matrix(get_xyz(mut_a), 3) + c(5, -3, 2))

  q <- genm_superpose_prot(mut_a, target = moved)
  expect_lt(max(abs(get_xyz(q) - moved)), 1e-10)

  # The Hessian and the covariance must be the original ones, rotated: with the
  # coordinates ordered x1 y1 z1 x2 ..., the rotation of the whole protein is
  # block-diagonal, one copy of rot per node.
  rot_all <- kronecker(diag(nsites), rot)
  k <- get_kmat(mut_a)
  expect_lt(max(abs(get_kmat(q) - rot_all %*% k %*% t(rot_all))), 1e-9 * max(abs(k)))
  cmat <- get_cmat(mut_a)
  expect_lt(max(abs(get_cmat(q) - rot_all %*% cmat %*% t(rot_all))), 1e-9 * max(abs(cmat)))
  expect_equal(get_evalue(q), get_evalue(mut_a))

  # what does not depend on orientation is untouched
  expect_identical(q$enm, mut_a$enm)
  expect_identical(genm_v_min(q), genm_v_min(mut_a))

  # a protein without modes stays without modes
  no_modes <- prot_from_enm(mut_a_enm, get_xyz(wt))
  expect_identical(genm_superpose_prot(no_modes, target = moved)$nma, NA)

  expect_error(genm_superpose_prot(mut_a, target = moved[-1]), "length")
})

test_that("the existing analysis getters work on a genm_prot", {
  expect_equal(length(get_msf_site(mut_a)), nsites)
  expect_true(all(is.finite(get_msf_site(mut_a)) & get_msf_site(mut_a) > 0))
  # cmat is the pseudo-inverse of kmat
  k <- get_kmat(mut_a)
  expect_lt(max(abs(k %*% get_cmat(mut_a) %*% k - k)), 1e-8 * max(abs(k)))
  # the mutation changes the spectrum, unlike lfenm
  expect_gt(abs(ddg_tds(wt, mut_a)), 0)
})

test_that("prot_from_enm leaves no modes; genm_add_nma adds them and changes nothing else", {
  p <- prot_from_enm(mut_a_enm, get_xyz(wt))
  expect_identical(p$nma, NA)
  expect_error(get_evalue(p))
  expect_error(get_cmat(p))

  q <- genm_add_nma(p)
  expect_equal(get_nmodes(q), 3 * nsites - 6)
  expect_identical(q[names(q) != "nma"], p[names(p) != "nma"])
  expect_identical(class(q), class(p))
  # recomputing on a protein that has modes gives the same modes
  expect_identical(genm_add_nma(q), q)

  expect_error(genm_add_nma(mut_a_enm))
})
