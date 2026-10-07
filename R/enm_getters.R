# Functions to query prot object -------------------------------------------------------


#' Get various properties of a prot object
#'
#' These accessors are the supported way to read a `prot` object built by
#' [set_enm()]. They fall into three groups: the model parameters (`get_enm_param`),
#' the nodes and their geometry (`get_nsites`, `get_site`, `get_pdb_site`,
#' `get_bfactor`, `get_xyz`), and the network and its normal modes (`get_kmat`,
#' `get_cmat`, `get_mode`, `get_evalue`, `get_umat`, `get_nmodes`).
#'
#' Prefer them over reaching into the list directly: a `prot` also carries
#' derived components that `penm` maintains for speed, whose internal consistency
#' is the package's business, and whose layout may change. The accessors are the
#' interface; the list structure is not.
#'
#' @param prot A protein with its associated ENM model, obtained using `set_enm`
#'
#' @seealso [set_enm()], which builds the `prot` these read.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' get_nsites(wt)                  # number of ENM nodes
#' get_nmodes(wt)                  # 3 * nsites - 6 for a generic structure
#' get_enm_param(wt)               # the arguments set_enm() was called with
#'
#' # site indexes are sequential; pdb_site is the numbering in the PDB file
#' head(get_site(wt))
#' head(get_pdb_site(wt))
#'
#' # xyz is a flat vector of length 3 * nsites, not a matrix
#' length(get_xyz(wt)) == 3 * get_nsites(wt)
#'
#' # kmat is 3N x 3N; umat is 3N x nmodes
#' dim(get_kmat(wt))
#' dim(get_umat(wt))
#'
#' @name get_prot_property
#'
NULL






#' @rdname get_prot_property
#'
#' @return \code{get_enm_param}: list of ENM parameters
#'
#' @export
#'
get_enm_param <- function(prot) prot$param



#' Get ENM node type
#'
#' @param prot is a prot object
#' @return ENM node type ("ca" or "sc")
#'
#' @noRd
#'
get_enm_node <- function(prot)  prot$param$node


#' Get ENM model type
#'
#' @param prot is a prot object
#' @return ENM model type ("anm", "hnm0", "hnm", "ming_wall", "pfanm")
#'
#'
#' @noRd
#'
get_enm_model <- function(prot)  prot$param$model


#' Get ENM contact distance cut-off
#'
#' @param prot is a prot object
#' @return d_max, the cut-off used to build the ENM network
#'
#'
#' @noRd
#'
get_d_max <- function(prot) prot$param$d_max


#' @rdname get_prot_property
#'
#'
#' @return \code{get_nsites}: number of sites (ENM nodes)
#'
#' @export
#'
get_nsites <- function(prot) prot$nodes$nsites


#' @rdname get_prot_property
#'
#'
#' @return \code{get_site}: site indexes, from `1` to `nsites`
#'
#' @export
#'
get_site  <- function(prot) prot$nodes$site


#' @rdname get_prot_property
#'
#'
#' @return \code{get_pdb_site}: site indexes, pdb numeration (resno)
#'
#' @export
#'
get_pdb_site <- function(prot) prot$nodes$pdb_site

#' @rdname get_prot_property
#'
#'
#' @return \code{get_bfactor}: B-factors of the X-ray file for the ENM nodes
#'
#' @export
#'
get_bfactor <- function(prot) prot$nodes$bfactor

#' @rdname get_prot_property
#'
#'
#' @return \code{get_xyz}: A 3*N vector of xyz coordinates of the N ENM nodes
#'
#' @export
#'

get_xyz <- function(prot)  prot$nodes$xyz


#' Get ENM graph
#'
#' @param prot is a prot object
#' @return graph, a tibble containing the graph representation of the ENM
#'
#'
#' @noRd
#'
get_graph <- function(prot) prot$graph

#' Get the unit vectors eij
#'
#' @param prot is a prot object
#' @return a matrix of size nedges x 3, containing unit vectors eij for all edges
#'
#'
#' @noRd
#'
get_eij <- function(prot) prot$eij


#' @rdname get_prot_property
#'
#'
#' @return \code{get_kmat}: the \code{3N x 3N} ENM network (stiffness) matrix, N = nsites
#'
#' @export
#'
get_kmat <- function(prot) prot$kmat

#' @rdname get_prot_property
#'
#'
#' @return \code{get_mode}: a vector of normal-mode indexes
#'
#' @export
#'
get_mode <- function(prot) prot$nma$mode

#' @rdname get_prot_property
#'
#'
#' @return \code{get_evalue}: a vector of normal-mode eigenvalues
#'
#' @export
#'
get_evalue <- function(prot) prot$nma$evalue

#' @rdname get_prot_property
#'
#'
#' @return \code{get_umat}: the matrix of eigenvectors U, of size \code{3 nsites x nmodes}
#'
#' @export
#'
get_umat <- function(prot) prot$nma$umat

#' @rdname get_prot_property
#'
#'
#' @return \code{get_cmat}: the ENM covariance matrix, of size \code{3 nsites x 3 nsites}
#'
#' @export
#'
get_cmat <- function(prot) prot$nma$cmat


#' @rdname get_prot_property
#'
#'
#' @return \code{get_nmodes}: the number of normal modes
#'
#' @export
#'
get_nmodes <- function(prot) length(prot$nma$mode)
