#' Depreciated, use sdmrs instead
#'
#' @keywords internal
#' @export
#'
dmrs_simulation <- function(...) {
  stop("dmrs_simulation demoted, use sdmrs instead; note sdmrs options are max_max and mean_max")
}

#' Depreciated, use admrs instead
#'
#' @keywords internal
#' @export
#'
dmrs_analytical <- function(...) {
  stop("dmrs_analytical demoted, use admrs instead; note admrs options are max_max and mean_max")
}
