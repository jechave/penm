# Run mutational scan

# Load libraries and functions

library(bio3d)
library(tidyverse) # includes dplyr and ggplot2
source("R/asm_functions.R")
source("R/set_parameters.R")


# Load data --------------------------------------------------------------

# file paths and names
dataset_path <- "data"
pdb_path <- "data/repaired_pdbs" # pdb files repaired using FoldX (by Elisha)
pdb_prefix <- "RepairPDB_"
run_path <- "runs/mut_scan/run_1" # directory for this run

# dataset
dataset <- read_csv(file.path(dataset_path,"monomers_by_pdb.csv.gz")) %>%
  dplyr::select(pdb, chain, n_sites, pdb_site_active) %>%
  mutate(pdb_site_active =  strsplit(pdb_site_active," ")) %>%
  mutate(pdb_site_active = map(pdb_site_active, as.integer)) %>%
  slice(1:10)

# parameters
run_par = tibble(
  nmut_per_site = 20,
  v0 = 0, # in kcal/mol
  mut_sd_min = 2,
  mut_dl_sigma = 1,
  fit_dg_thr = 0,
  fix_model = "moran",
  fix_n_eff = 1,
  fix_mut_rate = 1,
  beta = beta_boltzmann()
)
param <- do.call(set_param, as.list(run_par[1,]))

# Mutational scan

nprot <- nrow(dataset)
for (p in seq(nprot)) { # loop over proteins
  dp <- dataset[p,]
  pdb <- dp$pdb
  chain <- dp$chain
  pdb_site_active <- dp$pdb_site_active[[1]]
  print(paste("pdb:",pdb,"nsites:", dp$n_sites))

  # read structure
  wt <- read_pdb_ca(pdb, chain, path = pdb_path, prefix = pdb_prefix)

  # set up wt
  wt <- add_site_indexes(wt, pdb_site_active) # add indexes
  wt <- add_enm(wt, param$enm) # add enm info to wt
  wt <- add_v0(wt, v0 = param$enm$v0) # TODO: this is a hack
  wt$energy <- energy(wt, ideal = wt)
  wt$fitness <- fitness(wt,param$fit)
  wt$log_fit <- log_fitness_prot(wt, param$fit)

  nmut_per_site <- run_par$nmut_per_site
  nmut <- wt$nsites * nmut_per_site
  mutant_list <- vector("list",length = nmut)

  m <- 0
  for (site_mut in wt$site) {
    for (mut_at_site in seq(nmut_per_site)) {
      m <- m + 1
      mut <- get_mutant_site(wt, site_mut, wt, wt, param)
      # structural measures:
      msd <- 3*mean((mut$xyz - wt$xyz)^2, na.rm = T)
      # save
      mutant_list[[m + 1]] <-  tibble(pdb,
                                      chain,
                                      mutation = m,
                                      site = site_mut,
                                      pdb_site = mut$pdb_site[[site_mut]],
                                      v_min_wt = wt$energy$v_min,
                                      dv_activation_wt = wt$energy$dv_activation,
                                      v_min = mut$energy$v_min,
                                      dv_activation = mut$energy$dv_activation,
                                      msd)
    }
  }
  mutant_tbl <- bind_rows(mutant_list)

  # extra structural info
  distance_ca <- distance_to_active(wt$xyz, wt$site_active)
  cn_ca <- cn_graph(wt$enm$graph, wt$nsites)
  str_tbl <- tibble(site = wt$site,distance_ca,cn_ca)

  # join ouput
  out_tbl <- full_join(mutant_tbl,str_tbl)

  # save output
  out_file <- file.path(run_path,paste0(pdb, "_", chain, ".csv.gz"))
  write.csv(out_tbl, file = gzfile(out_file),row.names = FALSE)
}
