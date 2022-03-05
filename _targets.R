## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

use_files = dir("data/data_input", pattern = "mzML", full.names = TRUE)
excel_files = dir("data/data_input", pattern = ".xlsx", full.names = TRUE)

pkg_sha = as.character(utils::packageDescription("FTMS.peakCharacterization")$GithubSHA1)

pkg_tar = tar_target(pkg, pkg_sha)

data_mzml = expand_grid(
  data_function = rlang::syms("reading_scans_tile_windows"),
  mzml = use_files
) %>%
  dplyr::mutate(name = dplyr::case_when(
    grepl("97Cpos", mzml) ~ "97lipid",
    grepl("49Cpos", mzml) ~ "49lipid",
    grepl("1_ECF", mzml) ~ "1ecf",
    grepl("2_ECF", mzml) ~ "2ecf"
  ))

method_data = expand_grid(
  method_function = (c("noperc_nonorm",
                                  "perc99_nonorm",
                                  "singlenorm",
                                  "intsinglenorm",
                                  "doublenorm",
                                  "filtersd")),
  data_names = paste0("data_", data_mzml$name)) %>%
  dplyr::mutate(name = gsub("data_", "", data_names),
                data_sym = rlang::syms(data_names),
                method_sym = rlang::syms(method_function))

data_tar = tar_map(
  values = data_mzml,
  names = "name",
  tar_target(data, data_function(mzml, pkg))
)
method_tar = tar_map(
  values = method_data,
  unlist = FALSE,
  names = c("method_function", "name"),
  tar_target(method, method_sym(data_sym)),
  tar_target(rsd, single_rsd(method)),
  tar_target(zip, write_peaks_for_assignment(method)),
  tar_target(assign, assign_files(zip))
)

rsd_tar = tar_combine(
  rsd_combine,
  method_tar[[2]],
  command = merge_rsd(!!!.x),
  iteration = "list"
)

# rsd_tar = tar_combine(
#   rsd_data,
#   values = ends_with("97lipid")
# )

one_off_tar = tar_plan(
  tar_target(frequency_conversion, plot_frequency_conversion(data_97lipid)),
  tar_target(peak_ordering, plot_peak_ordering(data_97lipid)),
  tar_target(sliding_regions, plot_sliding_window_density(data_97lipid)),
  tar_target(peak_fit_plot, plot_peak_fitting(data_97lipid)),

  tar_render(manuscript, "doc/peakcharacterization_manuscript.Rmd")
)

list(pkg_tar, data_tar, method_tar, rsd_tar, one_off_tar)


# start target method_noperc_nonorm_161212_unlabeledAAs_2_ECF
# error target method_noperc_nonorm_161212_unlabeledAAs_2_ECF
# trying to get slot "elementMetadata" from an object of a basic class ("list") with no slots
