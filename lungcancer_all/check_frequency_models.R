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

run_mzml_list_FTMS = function(mzml_files, json_files = NULL, progress = TRUE, save_loc = ".", ...){

  mzml_exist = file.exists(mzml_files)

  if (!any(mzml_exist)) {
    stop("At least one of the mzML files supplied doesn't exist!")
  }

  if (!is.null(json_files)) {
    json_exist = file.exists(json_files)
    if (!any(json_exist)) {
      stop("At least one of the json files supplied doesn't exist!")
    }

    mzml_strip = gsub("mzML$", "", basename(mzml_files))
    json_strip = gsub("json$", "", basename(json_files))

    if (!all(all.equal(mzml_strip, json_strip))) {
      warning("Some of the mzML and json files don't match, are you sure they are all correct?")
    }
  }

  tictoc::tic()
  n_files = length(mzml_files)
  zip_files = purrr::map_chr(seq(1, n_files), function(i_file){

    in_mzml = mzml_files[i_file]

    if (!is.null(json_files)) {
      in_json = json_files[i_file]
    } else {
      in_json = NULL
    }

    zip_file = file.path(save_loc, gsub("mzML$", "zip", basename(in_mzml)))

    if (progress) {
      message(basename(in_mzml))
    }

    if (!file.exists(zip_file)) {
      FTMS.peakCharacterization:::log_message(paste0("Sample: ", basename(zip_file)))
      char_ms = CharacterizeMSPeaks$new(in_mzml, metadata_file = in_json, out_file = zip_file, ...)
      result = try({char_ms$run_all()})
    } else {
      print("file alreay exists!")
    }
    if (class(result) %in% "try-error") {
      out_result = result
    } else {
      out_result = zip_file
    }
    out_result
  })
  names(zip_files) = mzml_files
  tmp = tictoc::toc()
  saveRDS(zip_files, file = file.path(save_loc, "mzml_files_processed.rds"))

  message(paste0("processed: ", sum(!grepl("Error", zip_files))))
  message(paste0("errors: ", sum(grepl("Error", zip_files))))
  list(time = tmp, zip_files = zip_files)
}


peak_pick_samples <- read.table(here::here("lungcancer_all", "file_sample_info.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
all_mzml_files <- dir("/mlab/data/archives/FTMS_raw/CESB/mzml_data", recursive = TRUE, full.names = TRUE, pattern = "mzML$")

all_mzml_files <- data.frame(mzml_file = basename(all_mzml_files), mzml_path = all_mzml_files, stringsAsFactors = FALSE, sample = gsub(".mzML", "", basename(all_mzml_files)))

peak_pick_samples <- dplyr::left_join(peak_pick_samples, all_mzml_files, by = "sample")

files_04 = dir("/mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2022-04-02", pattern = "zip$")

files_05 = dir("/mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2022-05-11", pattern = "zip$")

missing_samples = gsub(".zip$", "", setdiff(files_04, files_05))

peak_pick_missing = dplyr::filter(peak_pick_samples, sample %in% missing_samples)

use_files <- unique(peak_pick_missing$mzml_path)

library(furrr)
plan(multicore)
FTMS.peakCharacterization::set_internal_map(furrr::future_map)

json_files = gsub(".mzML", ".json", use_files)

zip_ftms_dir = file.path("/mlab/scratch/cesb_data/zip_files/check_frequency_models", "FTMS_old_model")
dir.create(zip_ftms_dir)

frequency_model = c(0, -1/2, -1/3)
peak_region_finder = PeakRegionFinder$new(frequency_fit_description = frequency_model)
zip_ftms = run_mzml_list_FTMS(use_files, json_files = json_files, save_loc = zip_ftms_dir, raw_scan_filter = scan_time_rtime_filter, peak_finder = peak_region_finder)

saveRDS(zip_ftms, file = file.path(zip_ftms_dir, "zip_data.rds"))
textme::textme("old done")

library(ScanCentricPeakCharacterization)
# using the new model, in the new package to see if its
# the package or the model in general
run_mzml_list_SC = function(mzml_files, json_files = NULL, progress = TRUE, save_loc = ".", ...){

  mzml_exist = file.exists(mzml_files)

  if (!any(mzml_exist)) {
    stop("At least one of the mzML files supplied doesn't exist!")
  }

  if (!is.null(json_files)) {
    json_exist = file.exists(json_files)
    if (!any(json_exist)) {
      stop("At least one of the json files supplied doesn't exist!")
    }

    mzml_strip = gsub("mzML$", "", basename(mzml_files))
    json_strip = gsub("json$", "", basename(json_files))

    if (!all(all.equal(mzml_strip, json_strip))) {
      warning("Some of the mzML and json files don't match, are you sure they are all correct?")
    }
  }

  tictoc::tic()
  n_files = length(mzml_files)
  zip_files = purrr::map_chr(seq(1, n_files), function(i_file){

    in_mzml = mzml_files[i_file]

    if (!is.null(json_files)) {
      in_json = json_files[i_file]
    } else {
      in_json = NULL
    }

    zip_file = file.path(save_loc, gsub("mzML$", "zip", basename(in_mzml)))

    if (progress) {
      message(basename(in_mzml))
    }

    if (!file.exists(zip_file)) {
      ScanCentricPeakCharacterization:::log_message(paste0("Sample: ", basename(zip_file)))
      char_ms = SCCharacterizePeaks$new(in_mzml, metadata_file = in_json, out_file = zip_file, ...)
      result = try({char_ms$run_all()})
    } else {
      print("file alreay exists!")
    }
    if (class(result) %in% "try-error") {
      out_result = result
    } else {
      out_result = zip_file
    }
    out_result
  })
  names(zip_files) = mzml_files
  tmp = tictoc::toc()
  saveRDS(zip_files, file = file.path(save_loc, "mzml_files_processed.rds"))

  message(paste0("processed: ", sum(!grepl("Error", zip_files))))
  message(paste0("errors: ", sum(grepl("Error", zip_files))))
  list(time = tmp, zip_files = zip_files)
}
ScanCentricPeakCharacterization::set_internal_map(furrr::future_map)
zip_sc_dir = file.path("/mlab/scratch/cesb_data/zip_files/check_frequency_models", "SC_new_model")
dir.create(zip_sc_dir)
zip_sc = run_mzml_list_SC(use_files, json_files = json_files, save_loc = zip_sc_dir)
saveRDS(zip_sc, file = file.path(zip_sc_dir, "zip_data.rds"))
textme::textme("new done")
