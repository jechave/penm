#' Old energy function
#'
#' Demoted function, use enm_v_min and enm_g_entropy separately instead
#'
energy <- function(...) {
  stop("function energy deleted; call enm_v_min, enm_g_entropy, etc. separately")
}

#' Depreciated, use mutate_enm instead
#'
enm_update <- function(...) {
  stop("enm_update has been renamed to mutate_enm")
}


#' Depreciated, use sdmrs instead
#'
dmrs_simulation <- function(...) {
  stop("dmrs_simulation demoted, use sdmrs instead; note sdmrs options are max_max and mean_max")
}

#' Depreciated, use admrs instead
#'
dmrs_simulation <- function(...) {
  stop("dmrs_analytical demoted, use admrs instead; note admrs options are max_max and mean_max")
}
