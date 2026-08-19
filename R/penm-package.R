#'@details
#'  The \code{penm} package includes functions to calculate various Elastic Network Models
#'     for proteins and perform normal mode analysis (\code{\link{set_enm}}), to obtain
#'     mutant proteins and the corresponding mutant ENMs by perturbing the wild type
#'     (\code{\link{get_mutant_site}}), and to measure the resulting differences between
#'     wild type and mutant in energy (\code{\link{delta_energy}}), structure
#'     (\code{\link{delta_structure_by_site}}, \code{\link{delta_structure_by_mode}}), and
#'     motion (\code{\link{delta_motion_by_site}}, \code{\link{delta_motion_by_mode}}).
#'
#'  Mutants are identified by \code{(ensemble, site_mut, mutation)}; see
#'     \code{\link{penm_ensemble}} for what \code{ensemble} means and when it
#'     may be changed.
#'
#' @keywords internal
"_PACKAGE"

