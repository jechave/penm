# Join two independent mutational scans


# Libraries

library(plyr)
library(tidyverse)
library(bio3d)
library(fitdistrplus) # package to analyse and fit distributions
source("R/utility_functions.R")
source("R/analysis_functions.R")



# directory with mutant-scan files

path_1 <- "runs/mut_scan/run_4"
path_2 <- "runs/mut_scan/run_5"
path_output <- "runs/mut_scan/run_4_5"

# read data

files_1 <- dir(path = path_1, pattern = "csv.gz", full.names = TRUE, recursive = TRUE)
files_2 <- dir(path = path_2, pattern = "csv.gz", full.names = TRUE, recursive = TRUE)
if (length(files_1) != length(files_2)) {
  print("error: files1 != files2")
} else {
  mut_tab_1 <- tibble(file = files_1)
  mut_tab_1 <- mut_tab_1 %>%
    mutate(data = map(file, read.csv, as.is = TRUE, header = TRUE)) %>%
    unnest()

  mut_tab_2 <- tibble(file = files_2)
  mut_tab_2 <- mut_tab_2 %>%
    mutate(data = map(file, read.csv, as.is = TRUE, header = TRUE)) %>%
    unnest()
}


# remove the wt state from mut_tab_2

dat2 <- mut_tab_2 %>%
  filter(mutation != 0)

# change state identifier

shift <- max(mut_tab_1$mutation)
dat2 <- dat2 %>%
  mutate(mutation = mutation + shift)

# fix file names

dat1 <- mut_tab_1 %>%
  mutate(file = gsub(pattern = path_1, replacement = path_output, x = file))

dat2 <- dat2 %>%
  mutate(file = gsub(pattern = path_2, replacement = path_output, x = file))


# join the two sets

mut_tab <- bind_rows(dat1, dat2) %>%
  arrange(pdb, chain, site, mutation)



# write output

wf <- function(data, file) {
  write.csv(data, file = gzfile(file), row.names = FALSE)
}

dat <- mut_tab %>%
  group_by(file) %>%
  nest()

map2(dat$data, dat$file, wf)



