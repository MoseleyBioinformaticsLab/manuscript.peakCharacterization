library(drake)
library(FTMS.peakCharacterization)
library(rprojroot)
library(furrr)
plan(multiprocess)
project_root <- find_root(is_rstudio_project)
source("data_analysis/R/processing_functions.R")

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



# ---- Run it all!
full_plan <- rbind(setup_plan,
                   splitting_region_plan,
                   normalization_plan,
                   write_for_assignments,
                   write_for_other_analyses)

make(full_plan)
