#' Validate an ensemble label
#'
#' Errors unless \code{ensemble} is a single, non-missing, integer-valued
#' number (\code{1024} typed at the console, a double, is accepted; \code{1.5}
#' is not).
#'
#' @param ensemble The value to check.
#'
#' @return \code{ensemble}, invisibly, or an error.
#'
#' @noRd
check_ensemble <- function(ensemble) {
  # ensemble names which realization of the mutational process a mutant belongs
  # to (see ?penm_ensemble), so a malformed value is not a small problem: it
  # selects a realization nobody chose, and does so reproducibly, which is what
  # makes it dangerous -- the mistake never surfaces.
  #
  # Every bad value used to get one, because the key is built with paste(),
  # which stringifies anything. Measured before this check was added: NULL
  # yielded the key "-80-1" and hashed fine, and NA, "banana" and c(1, 2) all
  # returned a plausible integer.
  #
  # The length() != 1L test must precede is.na(), so the latter is never handed
  # a vector. ensemble != trunc(ensemble) rejects 1.5 while accepting 1024 typed
  # at the console, which is a double.
  if (is.null(ensemble) || length(ensemble) != 1L || !is.numeric(ensemble) ||
      is.na(ensemble) || ensemble != trunc(ensemble)) {
    stop("`ensemble` must be a single non-missing integer. ",
         "It names which realization of the mutational process a mutant ",
         "belongs to; see ?penm_ensemble.", call. = FALSE)
  }
  invisible(ensemble)
}

#' Map a mutant's identity to an RNG seed
#'
#' A mutant is identified by the tuple \code{(ensemble, site_mut, mutation)}:
#' \code{(site_mut, mutation)} names a specific mutation at a specific site,
#' and \code{ensemble} says which realization of the mutational process those
#' names refer to (see \code{?penm_ensemble}). The tuple is hashed, and the
#' hash is the seed: \code{ensemble} is a label, not a seed.
#'
#' The seed depends on the tuple alone, so a given \code{(site_mut, mutation)}
#' gets the same seed in every scan or trajectory that refers to it, however
#' many other mutations that scan includes.
#'
#' @param ensemble An integer naming the realization (see \code{?penm_ensemble}).
#' @param site_mut The mutated site (sequential index, not pdb_site).
#' @param mutation The mutation index at that site.
#'
#' @return An integer seed suitable for \code{set.seed()}.
#'
#' @noRd
mut_seed <- function(ensemble, site_mut, mutation) {
  # This function is the only place where a seed in the set.seed() sense exists.
  #
  # The tuple is hashed rather than packed arithmetically. Arithmetic keys
  # collide structurally: the previous scheme seed + site_mut * mutation gave
  # every divisor pair of the same product an identical random stream (for 228
  # sites x 10 mutations, 1710 of 2280 mutants shared a seed, up to 8 on one
  # value), and a positional key seed + site_mut * K + mutation merely moves the
  # collision up a level, making seed + 1 indistinguishable from mutation + 1.
  # Hashing makes distinctness a property of the hash instead of a property of a
  # hand-checked formula.
  #
  # Because the key never mentions nmut, an nmut = 50 run is a strict superset
  # of the nmut = 10 run.
  #
  # There is deliberately no separate ensemble-index slot beside a scan label.
  # The key once carried both, because sdmrs needs two independent ensembles
  # and, under the arithmetic key, scaling the label (1*seed, 2*seed) gave sets
  # that overlapped whenever nsites * nmut > seed. Under the hash two different
  # labels are already disjoint over a whole scan -- measured over 228 sites x
  # 10 mutations, the overlap is 0 for 1024 vs 1025 and for 1024 vs 7 -- so a
  # second slot was a second name for one axis. Independent ensembles come from
  # two different values of ensemble.

  # Checked here as well as at the set_enm() boundary: a direct penm:::
  # caller reaches this function without passing through it.
  check_ensemble(ensemble)
  key <- paste(ensemble, site_mut, mutation, sep = "-")
  # Use all 32 hash bits, then fold into the signed range set.seed() accepts
  # (it errors at or above 2^31). Truncating the hex string instead would throw
  # away bits and raise the collision rate: 7 hex digits is only 28 bits, which
  # for 10,000 mutants collides ~17% of the time versus ~1% at 32 bits.
  hex <- digest::digest(key, algo = "xxhash32")
  # strtoi() returns NA on 8 hex digits (it overflows signed int), so read the
  # halves separately and combine in double precision before folding.
  h <- strtoi(substr(hex, 1L, 4L), 16L) * 65536 + strtoi(substr(hex, 5L, 8L), 16L)
  as.integer(h %% 2147483647)
}

#' Evaluate an expression under a fixed seed, leaving the caller's RNG alone
#'
#' Evaluates \code{expr} after \code{set.seed(seed)}, then restores the
#' caller's \code{.Random.seed}, or removes it if there was none. The values
#' drawn inside \code{expr} are those \code{set.seed(seed)} gives.
#'
#' @param seed An integer seed, as returned by \code{mut_seed()}.
#' @param expr An expression to evaluate.
#'
#' @return The value of \code{expr}.
#'
#' @noRd
with_mut_seed <- function(seed, expr) {
  # set.seed() writes .Random.seed in the global environment, so seeding a
  # mutant's perturbations would also silently reseed whatever the caller was
  # doing: a loop that draws a mutant and then draws something of its own would
  # get its stream reset on every iteration.
  #
  # Written by hand rather than with withr::with_seed() to avoid taking a
  # dependency for one call site.
  #
  # If .Random.seed does not exist yet (no RNG use in the session so far) it is
  # removed again afterwards, so the session is returned to the state it was
  # actually in. That branch is deliberately not covered by a test: arranging
  # "no .Random.seed exists" means deleting it from the global environment, and
  # testthat would carry that side effect into later test files.
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  # expr is a promise, forced here: inside the function, and before on.exit
  # fires. Do not "simplify" this with force() or eval().
  expr
}
