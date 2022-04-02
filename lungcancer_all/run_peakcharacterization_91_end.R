library(methods)
library(FTMS.peakCharacterization)
pkg = utils::packageDescription("FTMS.peakCharacterization")
# actually run this stuff

scan_time_rtime_filter = function(raw_ms, min_time_difference = 4, rtime_limit = 7.5*60){
  scan_times = raw_ms$ms_info
  scan_times = scan_times[scan_times$scan %in% raw_ms$scan_range, ]

  scan_times = scan_times %>%
    dplyr::filter(rtime < rtime_limit)

  scan_times <- dplyr::mutate(scan_times, lag = rtime - dplyr::lag(rtime), lead = dplyr::lead(rtime) - rtime)

  high_lag <- scan_times$lag >= min_time_difference
  high_lag[is.na(high_lag)] <- TRUE
  high_lead <- scan_times$lead >= min_time_difference
  high_lead[is.na(high_lead)] <- TRUE

  na_lead_high_lag <- is.na(scan_times$lead) & high_lag
  na_lag_high_lead <- is.na(scan_times$lag) & high_lead

  keep_scans <- (na_lead_high_lag | high_lag) & (na_lag_high_lead | high_lead)
  raw_ms$set_scans(scan_range = scan_times$scan[keep_scans])
  raw_ms
}


peak_pick_samples <- read.table(here::here("lungcancer_all", "file_sample_info.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
all_mzml_files <- dir("/mlab/data/archives/FTMS_raw/CESB/mzml_data", recursive = TRUE, full.names = TRUE, pattern = "mzML$")

all_mzml_files <- data.frame(mzml_file = basename(all_mzml_files), mzml_path = all_mzml_files, stringsAsFactors = FALSE, sample = gsub(".mzML", "", basename(all_mzml_files)))

peak_pick_samples <- dplyr::left_join(peak_pick_samples, all_mzml_files, by = "sample")

curr_date = Sys.Date()
use_date = "2022-03-08"
if (curr_date != use_date) {
  warning("current date and the use date dont match!")
}
zip_save <- file.path("/mlab/scratch/cesb_data/zip_files", paste0("lung_matched_tissue-", use_date))
if (!dir.exists(zip_save)) {
  dir.create(zip_save)
}

existing_files <- dir(zip_save, pattern = ".zip", full.names = TRUE)
existing_samples <- gsub(".zip", "", basename(existing_files))

if (!file.exists(file.path(zip_save, "pkg.rds"))) {
  saveRDS(pkg, file.path(zip_save, "pkg.rds"))
} else {
  old_pkg = readRDS(file.path(zip_save, "pkg.rds"))
  curr_sha = digest::digest(pkg, algo = "sha256")
  old_sha = digest::digest(old_pkg, algo = "sha256")
  if (curr_sha != old_sha) {
    warning("saved pkg and current pkg don't match!")
  }
}

peak_pick_samples <- peak_pick_samples[!(peak_pick_samples$sample %in% existing_samples), ]

n_sample <- floor(nrow(peak_pick_samples) / 2)
use_files <- peak_pick_samples$mzml_path[(n_sample + 1):(nrow(peak_pick_samples))]

library(furrr)
plan(multicore)
set_internal_map(furrr::future_map)

json_files = gsub(".mzML", ".json", use_files)
zip_results = run_mzml_list(use_files, json_files = json_files, save_loc = zip_save, raw_scan_filter = scan_time_rtime_filter)

saveRDS(zip_results, file = file.path(zip_save, "zip_91_end.rds"))
