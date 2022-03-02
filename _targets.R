## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

use_files = dir(here::here("data/data_input"), pattern = "mzML", full.names = TRUE)
excel_files = dir(here::here("data/data_input"), pattern = ".xlsx", full.names = TRUE)

method_mzml = expand_grid(
  method_function = rlang::syms(c("noperc_nonorm",
                            "perc99_nonorm",
                            "singlenorm",
                            "singlenorm_int",
                            "doublenorm",
                            "filtersd")),
  mzml = use_files)
method_mzml = method_mzml %>%
  dplyr::mutate(rest_id = gsub(".mzML", "", basename(mzml)))
## tar_plan supports drake-style targets and also tar_target()
pkg = tar_target(pkg, utils::packageDescription('FTMS.peakCharacterization'))
method_peaks = tar_map(
  unlist = FALSE,
  values = method_mzml,
  names = c("method_function", "rest_id"),
  tar_target(data, reading_scans_tile_windows(mzml, pkg)),
  tar_target(processing, method_function(data)),
  tar_target(zip, write_peaks_for_assignment(processing)),
  tar_target(assign, assign_files(zip))
)

combine_scan_height = tar_combine(
  scan_height_cor,
  method_peaks[["data"]],
  command = correlate_scan_height(!!!.x)
)

one_offs = tar_plan(
  tar_target(frequency_conversion, plot_frequency_conversion(data_filtersd_97Cpos)),
  tar_target(peak_ordering, plot_peak_ordering(data_filtersd_97Cpos)),
  tar_target(sliding_regions, plot_sliding_window_density(data_filtersd_97Cpos)),
  tar_target(peak_fit_plot, plot_peak_fitting(data_filtersd_97Cpos)),

  tar_render(manuscript, "doc/peakcharacterization_manuscript.Rmd")
)

list(pkg, method_peaks, combine_scan_height, one_offs)
