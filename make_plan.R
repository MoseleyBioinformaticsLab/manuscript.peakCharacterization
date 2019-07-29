library(FTMS.peakCharacterization)
library(furrr)
future::plan(multicore(workers = 8))
set_internal_map(furrr::future_map)
library(rprojroot)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
library(patchwork)
library(dplyr)
library(glue)
library(smirfeTools)
project_root <- find_root(is_rstudio_project)

options(future.globals.maxSize = 2000 * 1024^2)
source("data_analysis/R/processing_functions.R")
source("data_analysis/R/graph_generation.R")

library(drake)
peak_pkg_description = utils::packageDescription('FTMS.peakCharacterization')

use_files = dir("data_analysis/data_input", pattern = "mzML", full.names = TRUE)

write_loc = "data_analysis/data_output"

analysis_plan = drake_plan(
  # check the version of peak characterization we are using, so if it changes, we rerun
  pkg = peak_pkg_description,
  data = target(
    file_in(input) %>% reading_scans_tile_windows(., pkg),
    transform = map(input = !!use_files, .id = FALSE)
  ),
  group1_method = target(
    group1_characterization(data),
    transform = map(data)
  ),
  group2_method = target(
    group2_characterization(data),
    transform = map(data)
  ),
  group3_method = target(
    group3_characterization(data),
    transform = map(data)
  ),
  group4_method = target(
    group4_characterization(data),
    transform = map(data)
  ),
  group5_method = target(
    group5_characterization(data),
    transform = map(data)
  ),
  group6_method = target(
    group6_characterization(data),
    transform = map(data)
  ),
  write_groups = target(
    write_peaks_for_assignment(group6_method),
    transform = map(group6_method)
  ),
  assign_groups = target(
    assign_files(write_groups),
    transform = map(write_groups)
  ),
  raw_assignments = target(
    read_assignments(assign_groups),
    transform = map(assign_groups)
  ),
  interesting_peaks = target(
    find_interesting_peaks(raw_assignments),
    transform = map(raw_assignments)
  )
)

make(analysis_plan, lock_envir = FALSE, memory_strategy = "memory")
