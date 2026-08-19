#' Validate an ensemble label
#'
#' \code{ensemble} names which realization of the mutational process a mutant
#' belongs to (see \code{?penm_ensemble}), so a malformed value is not a small
#' problem: it selects a realization nobody chose, and does so reproducibly,
#' which is what makes it dangerous -- the mistake never surfaces.
#'
#' Every bad value used to get one, because the key is built with
#' \code{paste()}, which stringifies anything. Measured before this check was
#' added: \code{NULL} yielded the key \code{"-80-1"} and hashed fine, and
#' \code{NA}, \code{"banana"} and \code{c(1, 2)} all returned a plausible
#' integer.
#'
#' The \code{length() != 1L} test must precede \code{is.na()}, so the latter is
#' never handed a vector. \code{ensemble != trunc(ensemble)} rejects \code{1.5}
#' while accepting \code{1024} typed at the console, which is a double.
#'
#' @param ensemble The value to check.
#'
#' @return \code{ensemble}, invisibly, or an error.
#'
#' @noRd
check_ensemble <- function(ensemble) {
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
#' names refer to (see \code{?penm_ensemble}).
#'
#' This function is the only place where a seed in the \code{set.seed()} sense
#' exists. \code{ensemble} is not one: it is a label, hashed together with the
#' rest of the tuple, and it is the hash that seeds the RNG.
#'
#' The tuple is hashed rather than packed arithmetically. Arithmetic keys
#' collide structurally: the previous scheme \code{seed + site_mut * mutation}
#' gave every divisor pair of the same product an identical random stream (for
#' 228 sites x 10 mutations, 1710 of 2280 mutants shared a seed, up to 8 on one
#' value), and a positional key \code{seed + site_mut * K + mutation} merely
#' moves the collision up a level, making \code{seed + 1} indistinguishable
#' from \code{mutation + 1}. Hashing makes distinctness a property of the hash
#' instead of a property of a hand-checked formula.
#'
#' Because the key never mentions \code{nmut}, a given \code{(site, mutation)}
#' draws the same stream in every scan that refers to it: an \code{nmut = 50}
#' run is a strict superset of the \code{nmut = 10} run.
#'
#' There is deliberately no separate ensemble-index slot beside a scan label.
#' The key once carried both, because \code{sdmrs} needs two independent
#' ensembles and, under the arithmetic key, scaling the label (\code{1*seed},
#' \code{2*seed}) gave sets that overlapped whenever \code{nsites * nmut >
#' seed}. Under the hash two different labels are already disjoint over a whole
#' scan -- measured over 228 sites x 10 mutations, the overlap is 0 for 1024 vs
#' 1025 and for 1024 vs 7 -- so a second slot was a second name for one axis.
#' Independent ensembles come from two different values of \code{ensemble}.
#'
#' @param ensemble An integer naming the realization (see \code{?penm_ensemble}).
#' @param site_mut The mutated site (sequential index, not pdb_site).
#' @param mutation The mutation index at that site.
#'
#' @return An integer seed suitable for \code{set.seed()}.
#'
#' @noRd
mut_seed <- function(ensemble, site_mut, mutation) {
  # Checked here as well as at the get_mutant_site() boundary: penmscan and any
  # direct penm::: caller reach this function without passing through it.
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
#' \code{set.seed()} writes \code{.Random.seed} in the global environment, so
#' seeding a mutant's perturbations also silently reseeds whatever the caller
#' was doing. A loop that draws a mutant and then draws something of its own
#' gets its stream reset on every iteration.
#'
#' Restoring the previous \code{.Random.seed} on exit changes no drawn value
#' inside \code{expr} -- the seeding is identical, only the aftermath differs.
#' Written by hand rather than with \code{withr::with_seed()} to avoid taking a
#' dependency for one call site.
#'
#' \code{expr} is a promise, forced on the last line: inside the function, and
#' before \code{on.exit} fires. Do not "simplify" this with \code{force()} or
#' \code{eval()}.
#'
#' If \code{.Random.seed} does not exist yet (no RNG use in the session so far)
#' it is removed again afterwards, so the session is returned to the state it
#' was actually in. That branch is deliberately not covered by a test:
#' arranging "no \code{.Random.seed} exists" means deleting it from the global
#' environment, and testthat would carry that side effect into later test
#' files.
#'
#' @param seed An integer seed, as returned by \code{mut_seed()}.
#' @param expr An expression to evaluate.
#'
#' @return The value of \code{expr}.
#'
#' @noRd
with_mut_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  expr
}

#' Check that a set of mutant seeds contains no duplicates
#'
#' Hash distinctness is probabilistic, not structural: over 32 bits the chance
#' any two of 2280 mutants collide is about 0.06%, rising to roughly 1% at
#' 10,000 mutants and higher still beyond that. A collision would silently make
#' two mutants non-independent, which is the exact failure this scheme exists to
#' remove, so assert rather than assume.
#'
#' @param seeds An integer vector of seeds.
#'
#' @return \code{seeds}, invisibly, or an error.
#'
#' @noRd
check_seeds_distinct <- function(seeds) {
  dup <- anyDuplicated(seeds)
  if (dup > 0L) {
    stop("Seed collision: mutant seeds are not distinct (first repeat at index ",
         dup, "). Vary `ensemble` and re-run.", call. = FALSE)
  }
  invisible(seeds)
}
