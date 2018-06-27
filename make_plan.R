library(drake)
library(FTMS.peakCharacterization)
library(rprojroot)
library(furrr)
plan(multiprocess)
project_root <- find_root(is_rstudio_project)
source("data_analysis/R/processing_functions.R")

peak_pkg_description = utils::packageDescription('FTMS.peakCharacterization')

setup_100cpos <- drake_plan(
  mzml_100cpos = file_in('data_analysis/data_input/100Cpos/100Cpos.mzML'),
  tiled_100cpos = average_scans_tile_windows(mzml_100cpos, peak_pkg_description)
)

setup_ecf <- drake_plan(
  mzml_ecf = file_in('data_analysis/data_input/ecf_data/7_HsPlt_13Cglc_p5Uthrombin_160728_polar_ECF_2.mzML'),
  tiled_ecf = average_scans_tile_windows(mzml_ecf, peak_pkg_description)
)

full_plan <- rbind(setup_100cpos, setup_ecf)

make(full_plan)
