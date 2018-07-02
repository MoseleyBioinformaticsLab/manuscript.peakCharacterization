library(drake)
library(FTMS.peakCharacterization)
library(rprojroot)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
library(patchwork)
library(dplyr)
library(furrr)
plan(multiprocess)
project_root <- find_root(is_rstudio_project)
source("data_analysis/R/processing_functions.R")
source("data_analysis/R/graph_generation.R")

peak_pkg_description = utils::packageDescription('FTMS.peakCharacterization')

# Get the data in, filter the scans, and add the windowed regions ----------
setup_100cpos <- drake_plan(
  mzml_100cpos = file_in('data_analysis/data_input/100Cpos/100Cpos.mzML'),
  tiled_100cpos = average_scans_tile_windows(mzml_100cpos, peak_pkg_description)
)

setup_ecf <- drake_plan(
  mzml_ecf = file_in('data_analysis/data_input/ecf_data/7_HsPlt_13Cglc_p5Uthrombin_160728_polar_ECF_2.mzML'),
  tiled_ecf = average_scans_tile_windows(mzml_ecf, peak_pkg_description)
)

setup_plan <- rbind(setup_100cpos, setup_ecf)

# Reduce and split the window regions with and without noise ---------------
splitting_regions_hasnoise <- drake_plan(
  hasnoise_100cpos = split_with_noise(tiled_100cpos),
  hasnoise_ecf = split_with_noise(tiled_ecf)
)

splitting_regions_nonoise <- drake_plan(
  nonoise_100cpos = split_without_noise(tiled_100cpos),
  nonoise_ecf = split_without_noise(tiled_ecf)
)

splitting_region_plan <- rbind(splitting_regions_hasnoise,
                               splitting_regions_nonoise)

# characterize with / without normalization -----
donormalization <- drake_plan(
  hasnorm_hasnoise_100cpos = run_characterization(hasnoise_100cpos, TRUE),
  hasnorm_nonoise_100cpos = run_characterization(nonoise_100cpos, TRUE),
  hasnorm_hasnoise_ecf = run_characterization(hasnoise_ecf, TRUE),
  hasnorm_nonoise_ecf = run_characterization(nonoise_ecf, TRUE)
)

nonormlization <- drake_plan(
  nonorm_hasnoise_100cpos = run_characterization(hasnoise_100cpos),
  nonorm_nonoise_100cpos = run_characterization(nonoise_100cpos),
  nonorm_hasnoise_ecf = run_characterization(hasnoise_ecf),
  nonorm_nonoise_ecf = run_characterization(nonoise_ecf)
)

normalization_plan <- rbind(donormalization,
                            nonormlization)

# create peak output for assignments --------------------------------------
write_for_assignments <- drake_plan(
  peaks_hasnorm_hasnoise_100cpos = write_peaks_for_assignment(hasnorm_hasnoise_100cpos,
                                                              file_out("data_analysis/data_output/peaks_hasnorm_hasnoise_100cpos.json")),
  peaks_hasnorm_nonoise_100cpos = write_peaks_for_assignment(hasnorm_nonoise_100cpos,
                                                             file_out("data_analysis/data_output/peaks_hasnorm_nonoise_100cpos.json")),
  peaks_nonorm_hasnoise_100cpos = write_peaks_for_assignment(nonorm_hasnoise_100cpos,
                                                             file_out("data_analysis/data_output/peaks_nonorm_hasnoise_100cpos.json")),
  peaks_nonorm_nonoise_100cpos = write_peaks_for_assignment(nonorm_nonoise_100cpos,
                                                            file_out("data_analysis/data_output/peaks_nonorm_nonoise_100cpos.json")),
  peaks_hasnorm_hasnoise_ecf = write_peaks_for_assignment(hasnorm_hasnoise_ecf,
                                                              file_out("data_analysis/data_output/peaks_hasnorm_hasnoise_ecf.json")),
  peaks_hasnorm_nonoise_ecf = write_peaks_for_assignment(hasnorm_nonoise_ecf,
                                                             file_out("data_analysis/data_output/peaks_hasnorm_nonoise_ecf.json")),
  peaks_nonorm_hasnoise_ecf = write_peaks_for_assignment(nonorm_hasnoise_ecf,
                                                             file_out("data_analysis/data_output/peaks_nonorm_hasnoise_ecf.json")),
  peaks_nonorm_nonoise_ecf = write_peaks_for_assignment(nonorm_nonoise_ecf,
                                                            file_out("data_analysis/data_output/peaks_nonorm_nonoise_ecf.json")),
)

# save peak output for other analyses -------------------------------------
write_for_other_analyses <- drake_plan(
  rds_hasnorm_hasnoise_100cpos = saveRDS(hasnorm_hasnoise_100cpos,
                                                              file_out("data_analysis/data_output/peaks_hasnorm_hasnoise_100cpos.rds")),
  rds_hasnorm_nonoise_100cpos = saveRDS(hasnorm_nonoise_100cpos,
                                                             file_out("data_analysis/data_output/peaks_hasnorm_nonoise_100cpos.rds")),
  rds_nonorm_hasnoise_100cpos = saveRDS(nonorm_hasnoise_100cpos,
                                                             file_out("data_analysis/data_output/peaks_nonorm_hasnoise_100cpos.rds")),
  rds_nonorm_nonoise_100cpos = saveRDS(nonorm_nonoise_100cpos,
                                                            file_out("data_analysis/data_output/peaks_nonorm_nonoise_100cpos.rds")),
  rds_hasnorm_hasnoise_ecf = saveRDS(hasnorm_hasnoise_ecf,
                                                          file_out("data_analysis/data_output/peaks_hasnorm_hasnoise_ecf.rds")),
  rds_hasnorm_nonoise_ecf = saveRDS(hasnorm_nonoise_ecf,
                                                         file_out("data_analysis/data_output/peaks_hasnorm_nonoise_ecf.rds")),
  rds_nonorm_hasnoise_ecf = saveRDS(nonorm_hasnoise_ecf,
                                                         file_out("data_analysis/data_output/peaks_nonorm_hasnoise_ecf.rds")),
  rds_nonorm_nonoise_ecf = saveRDS(nonorm_nonoise_ecf,
                                                        file_out("data_analysis/data_output/peaks_nonorm_nonoise_ecf.rds")),
)

# graph_results -----------------------------------------------------------
graph_results <- drake_plan(
  # raw spectrum graph
  # location of characterized peak
  # set of bad nap ratio peaks
  # offset vs M/Z and fit
  #
  # show the sliding window counts and the 99th percentile cutoff to denote noise
  sliding_regions_count_histogram = create_regions_count_histogram(tiled_100cpos,
                                                                   file_out("data_analysis/figure_output/sliding_regions_count_histogram.png")),
  # show a region after reduction that still has multiple peaks in it
  multiple_peak_region = create_multiple_peak_figure(hasnorm_nonoise_100cpos,
                                                     tiled_100cpos,
                                                     file_out("data_analysis/figure_output/multiple_peak_region_tiles.png"))
)


# ---- Run it all!
full_plan <- rbind(setup_plan,
                   splitting_region_plan,
                   normalization_plan,
                   write_for_assignments,
                   write_for_other_analyses,
                   graph_results)

make(full_plan)
