library(drake)
library(SIRM.FTMS.peakCharacterization)
library(rprojroot)
project_root <- find_root(is_rstudio_project)
source(file.path(project_root, "data_analysis/R/processing_functions.R"))

peak_pkg_description = utils::packageDescription('SIRM.FTMS.peakCharacterization')

peak_picking_plan <- drake_plan(
  zip_file = file_in(file.path(project_root, 'data_analysis/data_input/100Cpos.zip')),
  zip_data = raw_peakpicking(zip_file, peak_pkg_description)
)

noise_and_offsets_plan <- drake_plan(
  keep_noise_no_offset = keep_noise_no_offset_correspondence(zip_data),
  remove_noise_no_offset = remove_noise_remove_noise_no_offset_correspondence(zip_data),
  remove_noise_do_offset = remove_noise_do_offset_correspondence(zip_data)
)

full_plan <- rbind(peak_picking_plan, noise_and_offsets_plan)

make(full_plan)
