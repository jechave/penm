#TODO: run step by step and check results at each step
#TODO: test
#TODO: consider not including stress of covalent structure in v_min

# Run evolutionary simulation

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
run_path <- "runs/evo_trj/run_2" # directory for this run

# dataset
dataset <- read_csv(file.path(dataset_path,"monomers_by_pdb.csv.gz")) %>%
  select(pdb, chain, n_sites, pdb_site_active) %>%
  mutate(pdb_site_active =  strsplit(pdb_site_active," ")) %>%
  mutate(pdb_site_active = map(pdb_site_active, as.integer)) %>%
  filter(pdb == first(pdb)) # pick first protein for test runs

# parameters
run_par = tibble(
  nmut_per_site = 30,
  v0 = -30, # in kcal/mol
  mut_sd_min = 4,
  mut_dl_sigma = .2,
  fit_dg_thr = -5,
  fix_model = "mc",
  fix_mut_rate = 1.0,
  fix_n_eff = 1.e-10,
  beta = beta_boltzmann()
)
param <- do.call(set_param, as.list(run_par[1,]))


# Run evolutionary trajectories -----------------------------------------------

nprot <- nrow(dataset)
for (p in seq(nprot)) { # loop over proteins
  dp <- dataset[p,]
  pdb <- dp$pdb
  chain <- dp$chain
  pdb_site_active <- dp$pdb_site_active[[1]]
  print(paste("pdb:",pdb,"nsites:", dp$n_sites))

  # read structure
  wt_0 <- read_pdb_ca(pdb, chain, path = pdb_path, prefix = pdb_prefix)

  # set up wt_0
  wt_0 <- add_site_indexes(wt_0, pdb_site_active) # add indexes
  wt_0 <- add_enm(wt_0, param$enm) # add enm info to wt_0
  wt_0 <- add_v0(wt_0, v0 = param$enm$v0) # TODO: this is a hack
  wt_0$energy <- energy(wt_0, ideal = wt_0)
  wt_0$fitness <- fitness(wt_0,param$fit)
  wt_0$log_fit <- log_fitness_prot(wt_0, param$fit) # TODO: merge with previous

  # initialize trajectory
  nmut <- wt_0$nsites * param$trj$nmut_per_site
  trj <- vector("list",length = (nmut + 1 ) )

  # initial state
  wt_t <- wt_0
  substitution <- 0

  trj[[1]] <- save_state_initial(wt_t) # set and save initial data

  for (m in seq(nmut)) {
    mutant <- get_mutant(wt_t, wt_t, param)
    mut <- mutant$mut
    site_mut <- mutant$site_mut

    # structural measures:
    dxyz2 <- 3*mean((mut$xyz - wt_t$xyz)^2, na.rm = T)
    msd <- 3*mean((mut$xyz - wt_0$xyz)^2, na.rm = T)

    # fixation
    mut_p_fix <- p_fix_prot(mut, wt_t, param$fix)
    accept <- runif(1) <= mut_p_fix

    if (accept) {
      # Only necessary to update structure in graph for accepted mutations
      mut$eij <- eij_edge(mut$xyz, mut$graph$i, mut$graph$j)
      wt_t <- mut
      substitution <- substitution + 1
    }

    # save

    trj[[m + 1]] <- save_state_mutant(mut)
  }

  trj_tbl <- bind_rows(trj)

  # extra structural info
  distance_ca <- distance_to_active(wt_0$xyz, wt_0$site_active)
  cn_ca <- cn_graph(wt_0$enm$graph, wt_0$nsites)
  str_tbl <- tibble(site = wt_0$site,distance_ca,cn_ca)

  # join ouput
  out_tbl <- full_join(trj_tbl,str_tbl)

  # save output
  out_file <- file.path(run_path,paste0(pdb, "_", chain, ".csv.gz"))
  write.csv(out_tbl, file = gzfile(out_file),row.names = FALSE)
}
