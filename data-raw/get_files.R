library(bio3d)
library(here)


## Fetch PDB files and split to chain A only PDB files
ids <- c("2acy_A")
get_pdb_files(ids)



get_pdb_files <- function(ids) {
  raw.files <- get.pdb(ids, path = here("data"))
  files <- pdbsplit(raw.files, ids, path = here("data"))
}
