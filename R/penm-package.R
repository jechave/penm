#' penm: build and perturb elastic network models of proteins
#'
#' Builds Elastic Network Models (ENMs) of proteins, perturbs them, and measures
#' how far the perturbed protein has moved from the original in energy, structure
#' and motion. The name means "perturb ENM": mutation is the perturbation
#' implemented here, but the design allows for others.
#'
#' @details
#'  The \code{penm} package includes functions to calculate various Elastic Network Models
#'     for proteins and perform normal mode analysis (\code{\link{set_enm}}), to obtain
#'     mutant proteins and the corresponding mutant ENMs by perturbing the wild type
#'     (\code{\link{get_mutant_site}}), and to measure the resulting differences between
#'     wild type and mutant in energy (\code{\link{delta_energy}}), structure
#'     (\code{\link{delta_structure_by_site}}, \code{\link{delta_structure_by_mode}}), and
#'     motion (\code{\link{delta_motion_by_site}}, \code{\link{delta_motion_by_mode}}).
#'
#'  How a protein mutates is chosen when it is built, with \code{mut_model} in
#'     \code{\link{set_enm}}: the linearly forced ENM (\code{"lfenm"}, the
#'     default), or a generalized ENM (\code{"genm"}) in which the force
#'     constants follow the rest lengths and a mutant's structure is the minimum
#'     of its energy. A genm mutant comes without normal modes;
#'     \code{\link{set_enm_nma}} adds them. \code{\link{superpose_prot}} puts
#'     a protein in the orientation of another structure, for comparing proteins
#'     site by site.
#'
#'  Mutants are identified by \code{(ensemble, site_mut, mutation)}; see
#'     \code{\link{penm_ensemble}} for what \code{ensemble} means and when it
#'     may be changed.
#'
#'  \code{wt} and \code{mut} are roles in a comparison, not two kinds of object: both
#'     are \code{prot} objects, and any \code{prot} may play either part. This is why
#'     the \code{delta_*} functions take both. An evolutionary trajectory exploits
#'     it, feeding each mutant back as the wild type of the next generation.
#'
#'  Start with \code{vignette("penm")}; \code{vignette("genm")} covers the
#'     generalized model.
#'
#' @references
#'  Echave J (2008). Evolutionary divergence of protein structure: the linearly
#'  forced elastic network model. \emph{Chemical Physics Letters} \strong{457}(4--6),
#'  413--416. \doi{10.1016/j.cplett.2008.04.042}
#'
#'  Echave J, Fernandez FM (2010). A perturbative view of protein structural
#'  variation. \emph{Proteins} \strong{78}(1), 173--180. \doi{10.1002/prot.22553}
#'
"_PACKAGE"
