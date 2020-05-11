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

