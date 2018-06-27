library(drake)
library(FTMS.peakCharacterization)
library(rprojroot)
library(furrr)
plan(multiprocess)
project_root <- find_root(is_rstudio_project)
source(file.path(project_root, "data_analysis/R/processing_functions.R"))

peak_pkg_description = utils::packageDescription('FTMS.peakCharacterization')

setup_100cpos <- drake_plan(
  mzml_100cpos = file_in(file.path(project_root, 'data_analysis/data_input/100Cpos/100Cpos.mzML')),
  tiled_100cpos = average_scans_tile_windows(mzml_100cpos, peak_pkg_description)
)

make(setup_100cpos)
