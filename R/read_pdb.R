#' Check pdb file issues
#'
#'  Checks if the pdb file exists and, if it does, some other issues
#'
#' @param pdb_id The pdb id (lower case)
#' @param chain pdb chain
#' @param path path of pdb files
#' @param prefix prefix of pdb file names
#'
#' @return If pdb file does not exist, returns NA, otherwise,
#'     it returns \code{list(missing_residues, n_sites, n_unique_resno)}
#'
#' @export
#'
#' @examples
#'
#' @family pdb input functions
#'
pdb_file_check <- function(pdb_id, chain = NA_character_, path, prefix = "") {
  pdb <- pdb_id
  pdb_path <- path
  pdb_file <- paste0(prefix, pdb, "_", chain, ".pdb")
  pdb_file <- file.path(pdb_path, pdb_file)
  if (!file.exists(pdb_file)) {
    print(paste("warning", pdb_file, "not found: returning NA"))
    return(NA)
  }

  # read structure
  pdb <- read.pdb(file = pdb_file) # input structure

  missing_residues = FALSE
  if (!inspect.connectivity(pdb)){
    missing_residues = TRUE
    print(paste("Warning: ",pdb_id,"has missing residues"))
  }

  n_sites <- pdb$atom %>%
    filter(elety == "CA") %>%
    nrow()

  n_unique_resno = length(unique(pdb$atom$resno))


  lst(missing_residues, n_sites, n_unique_resno )
}

#' Read pdb side-chains
#'
#' Reads pdb file, returns side-chain coordinates and bfactors
#'
#' @param pdb_id  protein's pdb id
#' @param chain protein chain
#' @param path path to pdb files
#' @param prefix prefix of pdb file names
#'
#' @return list of center-of-mass coordinates lst(xyz, pdb_site, bfactor)
#'
#' @export
#'
#' @examples
#'
#' @family pdb input functions
#'
read_pdb_sc <- function(pdb_id, chain = NA_character_, path, prefix = "") {
  pdb <- pdb_id
  pdb_path <- path
  pdb_file <- paste0(prefix, pdb, "_", chain, ".pdb")
  pdb_file <- file.path(pdb_path, pdb_file)
  if (!file.exists(pdb_file)) {
    print(paste("warning", pdb_file, "not found: returning NA"))
    return(NA)
  }

  # read structure
  pdb <- bio3d::read.pdb(file = pdb_file) # input structure

  if (!inspect.connectivity(pdb)){
    print(paste("Warning: ",pdb_id,"has missing residues"))
  }

  # calculate xyz coordinates of ca, cb, com, qm (qb by Micheletti), ql (qb by Levitt)
  r = residue.coordinates(pdb)
  b = residue.bfactors(pdb)

  pdb_site <- r$site

  xyz <-  r$com.xyz
  xyz_na <- is.na(xyz)
  xyz[xyz_na] <- r$m.xyz[xyz_na] #when center of mass is NA, use Michelettis approximation to Cbeta coordinates
  xyz_na <- is.na(xyz)
  xyz[xyz_na] <- r$ca.xyz[xyz_na] # use CA coordinates if neither com or micheletti are defined

  b.c <- b$b.c
  b.m <- b$b.m
  b.a <- b$b.a

  bfactor <- b.c # center-of-mass bfactors
  bfactor[is.na(bfactor)] <- b.m[is.na(bfactor)] # Use Micheletti's bfactor when com bfactor undefined'
  bfactor[is.na(bfactor)] <- b.a[is.na(bfactor)] # Use CA bfactor otherwise


  lst(xyz, pdb_site, bfactor)

}

#' Read pdb alpha carbons
#'
#' Reads pdb file, returns alpha-carbon coordinates and bfactors
#'
#' @param pdb_id  protein's pdb id
#' @param chain protein chain
#' @param path path to pdb files
#' @param prefix prefix of pdb file names
#'
#' @return list of center-of-mass coordinates lst(xyz, pdb_site, bfactor)
#'
#' @export
#'
#' @examples
#'
#' @family pdb input functions
#'
read_pdb_ca <- function(pdb_id, chain = NA_character_, path, prefix = "") {
  pdb <- pdb_id
  pdb_path <- path
  pdb_file <- paste0(prefix, pdb, "_", chain, ".pdb")
  pdb_file <- file.path(pdb_path, pdb_file)
  if (!file.exists(pdb_file)) {
    print(paste("warning", pdb_file, "not found: returning NA"))
    return(NA)
  }
  # read structure
  pdb <- read.pdb(file = pdb_file) # input structure
  # select CA
  sel <- atom.select(pdb, chain = chain, elety = "CA")
  # prepare output
  bfactor <- as.numeric(pdb$atom[sel$atom, c("b")])
  pdb_site <- as.numeric(pdb$atom[sel$atom, c("resno")])
  nsites <- length(pdb_site)
  xyz.calpha <-
    matrix(pdb$xyz[sel$xyz],
           ncol = nsites,
           nrow = 3,
           byrow = F)

  output <-
    list(
      "xyz" = as.vector(xyz.calpha),
      "pdb_site" = pdb_site,
      "bfactor" = bfactor
    )
  output
}
