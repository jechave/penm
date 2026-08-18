# Package data ------------------------------------------------------------

#' Acylphosphatase structure (PDB 2ACY, chain A)
#'
#' The crystal structure of bovine common-type acylphosphatase, chain A, as read by
#' [bio3d::read.pdb()]. Supplied so that the examples throughout `penm` can be run
#' without downloading anything.
#'
#' A small single-domain protein (98 residues), which keeps the normal-mode
#' calculations in the examples fast.
#'
#' @format A `pdb` object (class `c("pdb", "sse")`) as returned by
#'   [bio3d::read.pdb()], a list with components `atom`, `xyz`, `helix`, `sheet`,
#'   `calpha` and `call`. It holds 870 atoms over 98 residues, all in chain A.
#'
#' @source Protein Data Bank entry
#'   [2ACY](https://www.rcsb.org/structure/2ACY), chain A.
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, frustrated = FALSE)
#' get_nsites(wt)
"pdb_2acy_A"
