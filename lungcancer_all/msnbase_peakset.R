library(MSnbase)
source("./R/msnbase_xcms.R")
all_mzml_files <- dir("/mlab/data/archives/FTMS_raw/CESB/mzml_data", recursive = TRUE, full.names = TRUE, pattern = "mzML$")
names(all_mzml_files) = gsub(".mzML$", "", basename(all_mzml_files))

all_zip_files = dir("/mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2022-03-08", pattern = ".zip$", full.names = TRUE)
names(all_zip_files) = gsub(".zip$", "", basename(all_zip_files))

zip_mzml = intersect(names(all_mzml_files), names(all_zip_files))

all_mzml_files = all_mzml_files[zip_mzml]

msnbase_peaks = purrr::imap(all_mzml_files, function(in_mzml, id_mzml){
  msnbase_centroid(in_mzml, sample_id = id_mzml)
})

sample_id = gsub("_peaks.json", "", basename(all_json))
json_peaks = purrr::map(all_json, function(.x){
  tmp = jsonlite::fromJSON(.x)$Peaks
  sample = gsub("_peaks.json", "", basename(.x))
  tmp$sample = sample
  tmp
})

saveRDS(json_peaks, "data/data_output/lung_data/lung_xcalibur_peaks.rds")
