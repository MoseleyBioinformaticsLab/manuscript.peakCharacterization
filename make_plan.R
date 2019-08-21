library(FTMS.peakCharacterization)
library(furrr)
future::plan(multicore(workers = 8))
set_internal_map(furrr::future_map)
library(rprojroot)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
library(patchwork)
library(ggridges)
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
  # read in and prep each of the mzML files
  data = target(
    file_in(input) %>% reading_scans_tile_windows(., pkg),
    transform = map(input = !!use_files, .id = FALSE)
  ),
  # for each of the ways we could process the data, run it on each one of the
  # prepped data files
  method = target(
    method_fun(data),
    transform = cross(method_fun = c(noperc_nonorm,
                                     perc99_nonorm,
                                     singlenorm,
                                     singlenorm_int,
                                     doublenorm,
                                     filtersd), data)
  ),
  zip = target(
    write_peaks_for_assignment(method),
    transform = map(method)
  ),
  frequency_conversion = target(
    plot_frequency_conversion(data)
  ),
  peak_ordering = target(
    plot_peak_ordering(method_filtersd_data)
  ),
  sliding_regions = target(
    plot_sliding_window_density(data)
  ),
  # here we want to combine the ways a data file was processed and compare
  # their rsd's across processing methods.
  rsd = target(
    rsd_info(method),
    transform = combine(method, .by = data)
  ),
  rsd_plot = target(
    plot_rsd_differences(rsd),
    transform = map(rsd)
  ),
  rsd_table = target(
    summarize_rsd(rsd),
    transform = map(rsd)
  ),
  manuscript = rmarkdown::render(
    knitr_in("peakcharacterization_manuscript.Rmd"),
    output_file = file_out("peakcharacterization.docx"),
    quiet = TRUE
  )
)

make(analysis_plan, lock_envir = FALSE, memory_strategy = "preclean")
