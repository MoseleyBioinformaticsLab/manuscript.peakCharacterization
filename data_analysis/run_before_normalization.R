library(SIRM.FTMS.peakCharacterization)
library(rprojroot)
library(parallel)
options(mc.cores = 10)
project_root <- find_root(is_rstudio_project)
source(file.path(project_root, "data_analysis/R/processing_functions.R"))

peak_pkg_description = utils::packageDescription('SIRM.FTMS.peakCharacterization')

zip_files <- dir(file.path(project_root, "data_analysis", "data_input"), full.names = TRUE)

run_info <- mclapply(zip_files, function(in_file){
  zip_data <- raw_peakpicking(in_file, peak_pkg_description)
  zip_data$peak_finder$create_correspondent_peaks(median_corrected = FALSE)
  zip_data$peak_finder$collapse_correspondent_peaks()
  out_path <- file.path(project_root, "data_analysis", "data_output")
  out_file <- gsub(".zip", ".rds", basename(in_file))
  saveRDS(zip_data, file.path(out_path, out_file))
  file.path(out_path, out_file)
})
