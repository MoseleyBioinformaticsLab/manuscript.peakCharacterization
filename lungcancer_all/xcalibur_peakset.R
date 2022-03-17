library(jsonlite)
all_json = dir(here::here("ftms_artifacts/Peaklists"), full.names = TRUE)

sample_id = gsub("_peaks.json", "", basename(all_json))
json_peaks = purrr::map(all_json, function(.x){
  tmp = jsonlite::fromJSON(.x)$Peaks
  sample = gsub("_peaks.json", "", basename(.x))
  tmp$sample = sample
  tmp
})

saveRDS(json_peaks, "data/data_output/lung_data/lung_xcalibur_peaks.rds")
