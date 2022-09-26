## code to prepare `DATASET` dataset goes here

# usethis::use_data("DATASET")

library(tidyverse)
library(bio3d)
library(here)
library(jefuns)
library(penm)

load(here("data/pdb_2acy_A.rda"))

prot_2acy_A_ming_wall_ca <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE)

usethis::use_data(prot_2acy_A_ming_wall_ca, overwrite = TRUE)



