# Guards on mut_seed(): a mutant's identity must map to a distinct RNG stream.
#
# The scheme this replaced, seed + site_mut * mutation, is a product and so is
# not injective: every divisor pair of the same product shared a stream. These
# tests pin the properties that failure violated, plus the ones violated by the
# positional key (seed + site*K + mutation) considered as an alternative.

library(penm)

test_that("historically colliding pairs get distinct seeds and distinct draws", {
  # Under the old key these four all mapped to seed + 12 and drew the identical
  # rnorm() sequence, despite perturbing different edges.
  pairs <- list(c(2, 6), c(3, 4), c(4, 3), c(12, 1))

  seeds <- vapply(pairs, function(p) penm:::mut_seed(241956, p[1], p[2]), numeric(1))
  expect_equal(length(unique(seeds)), length(pairs))

  draws <- lapply(seeds, function(s) { set.seed(s); rnorm(5) })
  expect_equal(length(unique(draws)), length(pairs))
})

test_that("a full scan produces no colliding seeds", {
  # 228 sites x 10 mutations: the census that exposed the original bug, where
  # 2280 mutants collapsed onto 1134 distinct seeds.
  g <- expand.grid(j = seq(228), m = seq(10))
  keys <- mapply(function(j, m) penm:::mut_seed(241956, j, m), g$j, g$m)

  expect_equal(length(unique(keys)), nrow(g))
  expect_equal(anyDuplicated(keys), 0L)
})

test_that("seeds stay within the integer range set.seed() accepts", {
  # set.seed() errors above 2^31; a key that overflows would abort a scan.
  g <- expand.grid(j = seq(500), m = seq(50))
  keys <- mapply(function(j, m) penm:::mut_seed(241956, j, m), g$j, g$m)

  expect_true(all(keys <= .Machine$integer.max))
  expect_true(all(is.finite(keys)))
  expect_silent(set.seed(max(keys)))
})

test_that("a mutant's stream does not depend on nmut", {
  # (site, mutation) names a specific substitution, so scans that differ only
  # in how many mutants they draw must agree on the ones they share: an
  # nmut = 50 run is a strict superset of an nmut = 10 run.
  draw <- function() { set.seed(penm:::mut_seed(1024, 80, 3)); rnorm(3) }
  expect_equal(draw(), draw())

  keys10 <- vapply(seq(10), function(m) penm:::mut_seed(1024, 80, m), numeric(1))
  keys50 <- vapply(seq(50), function(m) penm:::mut_seed(1024, 80, m), numeric(1))
  expect_equal(keys10, keys50[seq(10)])
})

test_that("different seeds give disjoint scans", {
  # The regression that killed the positional key seed + site*K + mutation:
  # shifting the seed by 1 was indistinguishable from shifting mutation by 1,
  # so scans with seed = 1 and seed = 2 shared 90% of their streams.
  g <- expand.grid(j = seq(228), m = seq(10))
  k1 <- mapply(function(j, m) penm:::mut_seed(1, j, m), g$j, g$m)
  k2 <- mapply(function(j, m) penm:::mut_seed(2, j, m), g$j, g$m)

  expect_equal(length(intersect(k1, k2)), 0L)
})

test_that("the two sdmrs ensembles are disjoint", {
  # Previously separated by scaling the seed (1*seed, 2*seed), which overlapped
  # whenever nsites * nmut > seed.
  g <- expand.grid(j = seq(228), m = seq(10))
  e1 <- mapply(function(j, m) penm:::mut_seed(1024, j, m, 1L), g$j, g$m)
  e2 <- mapply(function(j, m) penm:::mut_seed(1024, j, m, 2L), g$j, g$m)

  expect_equal(length(intersect(e1, e2)), 0L)
})

test_that("check_seeds_distinct() fails loud on a collision", {
  expect_silent(penm:::check_seeds_distinct(c(1L, 2L, 3L)))
  expect_error(penm:::check_seeds_distinct(c(1L, 2L, 1L)), "Seed collision")
})
