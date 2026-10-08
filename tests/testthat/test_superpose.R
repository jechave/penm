# superpose_prot (R/superpose.R).
#
# Each protein is moved by a known rotation and translation, then superposed:
# what depends on orientation must come out rotated by exactly that rotation,
# and nothing else may change. Nothing here is frozen.

load(test_path("fixtures", "pdb_2acy_A.rda"))

# a rotation about an axis through the origin
rotation_about <- function(axis, angle) {
  u <- axis / sqrt(sum(axis^2))
  cross <- matrix(c(0, u[3], -u[2], -u[3], 0, u[1], u[2], -u[1], 0), 3)
  diag(3) * cos(angle) + sin(angle) * cross + (1 - cos(angle)) * tcrossprod(u)
}
rot <- rotation_about(c(1, 2, -1), 0.9)
shift <- c(5, -3, 2)
move <- function(xyz) as.vector(rot %*% matrix(xyz, nrow = 3) + shift)

wt_lfenm <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5)
mut_lfenm <- get_mutant_site(wt_lfenm, site_mut = 80, mutation = 3)
wt_genm <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
                   mut_model = "genm", d_max_graph = 14)
mut_genm <- set_enm_nma(get_mutant_site(wt_genm, site_mut = 80, mutation = 3))
proteins <- list(lfenm_wt = wt_lfenm, lfenm_mutant = mut_lfenm, genm_mutant = mut_genm)


test_that("superposing a moved structure onto itself recovers the rotation", {
  xyz <- get_xyz(wt_lfenm)
  moved <- as.vector(bio3d::fit.xyz(fixed = move(xyz), mobile = xyz,
                                    fixed.inds = seq_along(xyz), mobile.inds = seq_along(xyz)))
  r <- penm:::superposition_rotation(xyz, moved)
  expect_equal(r, rot)
  expect_equal(crossprod(r), diag(3))
  expect_equal(det(r), 1)
})

test_that("everything that depends on orientation is rotated, and nothing else changes", {
  for (name in names(proteins)) {
    p <- proteins[[name]]
    q <- superpose_prot(p, move(get_xyz(p)))
    rot_all <- kronecker(diag(get_nsites(p)), rot)

    expect_lt(max(abs(get_xyz(q) - move(get_xyz(p)))), 1e-9)
    k <- get_kmat(p)
    expect_lt(max(abs(get_kmat(q) - rot_all %*% k %*% t(rot_all))), 1e-9 * max(abs(k)))
    cmat <- get_cmat(p)
    expect_lt(max(abs(get_cmat(q) - rot_all %*% cmat %*% t(rot_all))), 1e-9 * max(abs(cmat)))
    # eigenvectors: rotated, each up to the sign convention, which still holds
    overlap <- abs(colSums(get_umat(q) * (rot_all %*% get_umat(p))))
    expect_lt(max(abs(overlap - 1)), 1e-9)
    expect_identical(penm:::canonical_sign(get_umat(q)), get_umat(q))
    if (!is.null(p$internal$eij)) {
      expect_lt(max(abs(q$internal$eij - p$internal$eij %*% t(rot))), 1e-12)
    }

    # independent of orientation: unchanged
    expect_identical(get_evalue(q), get_evalue(p))
    expect_identical(q$graph, p$graph)
    expect_identical(q$nodes[names(q$nodes) != "xyz"], p$nodes[names(p$nodes) != "xyz"])
    expect_identical(q$param, p$param)
    expect_equal(get_msf_site(q), get_msf_site(p))
    # dij is left as it was, and still the edge lengths of the moved structure
    expect_equal(penm:::dij_edge(get_xyz(q), q$graph$i, q$graph$j), q$graph$dij)
  }
})

test_that("superposing back onto the original structure restores the protein", {
  for (name in names(proteins)) {
    p <- proteins[[name]]
    back <- superpose_prot(superpose_prot(p, move(get_xyz(p))), get_xyz(p))
    expect_lt(max(abs(get_xyz(back) - get_xyz(p))), 1e-9)
    expect_lt(max(abs(get_kmat(back) - get_kmat(p))), 1e-9 * max(abs(get_kmat(p))))
    expect_lt(max(abs(get_cmat(back) - get_cmat(p))), 1e-9 * max(abs(get_cmat(p))))
  }
})

test_that("an lfenm protein mutates the same way in any orientation", {
  # moving the wild type, then mutating it, gives the mutant moved: this needs
  # kmat, the modes and eij all rotated consistently
  moved_wt <- superpose_prot(wt_lfenm, move(get_xyz(wt_lfenm)))
  mutant_of_moved <- get_mutant_site(moved_wt, site_mut = 80, mutation = 3)
  expect_gt(max(abs(get_xyz(mut_lfenm) - get_xyz(wt_lfenm))), 0.01)
  expect_lt(max(abs(get_xyz(mutant_of_moved) - move(get_xyz(mut_lfenm)))), 1e-9)
})

test_that("a protein without modes stays without modes, and the size is checked", {
  no_modes <- get_mutant_site(wt_genm, site_mut = 80, mutation = 3)
  expect_identical(superpose_prot(no_modes, move(get_xyz(no_modes)))$nma, NA)
  expect_error(superpose_prot(mut_genm, get_xyz(mut_genm)[-1]), "length")
})
