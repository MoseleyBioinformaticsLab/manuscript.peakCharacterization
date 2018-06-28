library(drake)
library(FTMS.peakCharacterization)
library(rprojroot)
library(furrr)
plan(multiprocess)
project_root <- find_root(is_rstudio_project)
source("data_analysis/R/processing_functions.R")

peak_pkg_description = utils::packageDescription('FTMS.peakCharacterization')

#----- Get the data in, filter the scans, and add the windowed regions
setup_100cpos <- drake_plan(
  mzml_100cpos = file_in('data_analysis/data_input/100Cpos/100Cpos.mzML'),
  tiled_100cpos = average_scans_tile_windows(mzml_100cpos, peak_pkg_description)
)

setup_ecf <- drake_plan(
  mzml_ecf = file_in('data_analysis/data_input/ecf_data/7_HsPlt_13Cglc_p5Uthrombin_160728_polar_ECF_2.mzML'),
  tiled_ecf = average_scans_tile_windows(mzml_ecf, peak_pkg_description)
)

setup_plan <- rbind(setup_100cpos, setup_ecf)

# ---- Reduce and split the window regions with and without noise
splitting_regions_hasnoise <- drake_plan(
  hasnoise_100cpos = split_with_noise(tiled_100cops),
  hasnoise_ecf = split_with_noise(tiled_ecf)
)

splitting_regions_nonoise <- drake_plan(
  nonoise_100cpos = split_without_noise(tiled_100cpos),
  nonoise_ecf = split_without_noise(tiled_ecf)
)

splitting_region_plan <- rbind(splitting_regions_hasnoise,
                               splitting_regions_nonoise)

# ---- characterize with and without normalization as well
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

# ---- Run it all!
full_plan <- rbind(setup_plan,
                   splitting_region_plan,
                   normalization_plan)

make(full_plan)
