#' Map a mutant's identity to an RNG seed
#'
#' A mutant is identified by the tuple \code{(seed, ensemble, site_mut,
#' mutation)}: \code{(site_mut, mutation)} names a specific substitution at a
#' specific site, \code{ensemble} distinguishes independent replicate ensembles
#' (used by \code{sdmrs}), and \code{seed} labels the whole scan.
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
#' @param seed An integer labelling the scan.
#' @param site_mut The mutated site (sequential index, not pdb_site).
#' @param mutation The mutation index at that site.
#' @param ensemble An integer distinguishing independent ensembles.
#'
#' @return An integer seed suitable for \code{set.seed()}.
#'
#' @noRd
mut_seed <- function(seed, site_mut, mutation, ensemble = 1L) {
  key <- paste(seed, ensemble, site_mut, mutation, sep = "-")
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
         dup, "). Vary `seed` and re-run.", call. = FALSE)
  }
  invisible(seeds)
}
