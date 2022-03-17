source("./packages.R")
lapply(list.files(here::here("./R"), full.names = TRUE), source)
library(knitrProgressBar)
all_mzml_files <- dir("/mlab/data/archives/FTMS_raw/CESB/mzml_data", recursive = TRUE, full.names = TRUE, pattern = "mzML$")
names(all_mzml_files) = gsub(".mzML$", "", basename(all_mzml_files))

all_zip_files = dir("/mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2022-03-08", pattern = ".zip$", full.names = TRUE)
names(all_zip_files) = gsub(".zip$", "", basename(all_zip_files))

zip_mzml = intersect(names(all_mzml_files), names(all_zip_files))

all_mzml_files = all_mzml_files[zip_mzml]

pb = knitrProgressBar::progress_estimated(length(all_mzml_files))
msnbase_peaks = purrr::imap(all_mzml_files, function(in_mzml, id_mzml){
  tmp = msnbase_centroid(in_mzml, sample_id = id_mzml)
  knitrProgressBar::update_progress(pb)
  tmp
})

saveRDS(msnbase_peaks, "data/data_output/lung_data/lung_msnbase_peaks.rds")
