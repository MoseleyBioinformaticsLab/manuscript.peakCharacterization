library(FTMS.peakCharacterization)
library(furrr)
future::plan(multicore)
set_internal_map(furrr::future_map)
library(rprojroot)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
library(patchwork)
library(dplyr)
project_root <- find_root(is_rstudio_project)
source("data_analysis/R/processing_functions.R")
source("data_analysis/R/graph_generation.R")

library(drake)
peak_pkg_description =

use_files = dir("data_analysis/data_input", pattern = "mzML", full.names = TRUE)

analysis_plan = drake_plan(
  # check the version of peak characterization we are using, so if it changes, we rerun
  pkg = target(
    # Triggers are always checked even though commands do not always run:
    trigger = trigger(change = utils::packageDescription('FTMS.peakCharacterization'))
  ),
  data = target(
    file_in(input) %>% reading_scans_tile_windows(., pkg),
    transform = map(input = !!use_files, .id = FALSE)
  ),
  final_method = target(
    final_characterization(data, pkg),
    transform = map(data)
  )
)

make(analysis_plan, lock_envir = FALSE, memory_strategy = "memory")
