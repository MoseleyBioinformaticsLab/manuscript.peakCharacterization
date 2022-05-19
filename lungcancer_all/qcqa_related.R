# this is where we are running some stuff to check where we are getting scan
# level QC/QA, besides in the actual peaks.
all_zip = dir("/mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2022-05-11", pattern = "zip$", full.names = TRUE)
library(FTMS.peakCharacterization)

scan_time_rtime_filter = function(scan_times, min_time_difference = 4, rtime_limit = 7.5*60){

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
  scan_times[keep_scans, ]
}

rawms_time_rtime_filter = function(raw_ms, min_time_difference = 4, rtime_limit = 7.5*60){
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


convert_scans2 = function(mz_df_list, frequency_fit_description, mz_fit_description){
  ... = NULL
  if (is.null(names(mz_df_list))) {
    names(mz_df_list) = purrr::map_chr(mz_df_list, ~ .x$scan[1])
  }
  mz_frequency = FTMS.peakCharacterization:::internal_map$map_function(mz_df_list, function(in_scan){
    #message(scan_name)
    out_scan = FTMS.peakCharacterization:::convert_mz_frequency(in_scan, ...)
    out_scan
  })
  #log_message("converted frequencies")
  frequency_fits = FTMS.peakCharacterization:::internal_map$map_function(mz_frequency, function(in_freq){
    use_peaks = in_freq$convertable
    tmp_fit = FTMS.peakCharacterization:::fit_exponentials(in_freq$mean_mz[use_peaks], in_freq$mean_frequency[use_peaks], frequency_fit_description)
    tmp_fit$scan = in_freq[1, "scan"]
    tmp_fit
  })
  #log_message("fit frequencies")

  mz_fits = FTMS.peakCharacterization:::internal_map$map_function(mz_frequency, function(in_freq){
    use_peaks = in_freq$convertable
    tmp_fit = FTMS.peakCharacterization:::fit_exponentials(in_freq$mean_frequency[use_peaks], in_freq$mean_mz[use_peaks], mz_fit_description)
    tmp_fit$scan = in_freq[1, "scan"]
    tmp_fit
  })
  #log_message("fit mz")

  frequency_coefficients = purrr::map_df(frequency_fits, function(.x){
    tmp_df = as.data.frame(matrix(.x$coefficients, nrow = 1))
    tmp_df$scan = .x$scan
    tmp_df
  })
  #log_message("extracted coefficients")

  mz_coefficients = purrr::map_df(mz_fits, function(.x){
    tmp_df = as.data.frame(matrix(.x$coefficients, nrow = 1))
    tmp_df$scan = .x$scan
    tmp_df
  })

  first_slope = min(which(frequency_fit_description != 0))
  bad_coefficients = boxplot.stats(frequency_coefficients[[first_slope]])$out

  return(list(coefficients = frequency_coefficients[[first_slope]], bad = bad_coefficients))

}

check_scans_removed = function(in_zip){
  zip_ms = ZipMS$new(in_zip)
  zip_ms$load_raw()
  zip_ms$load_peak_finder()

  raw_ms = zip_ms$raw_ms$clone(deep = TRUE)
  raw_ms = rawms_time_rtime_filter(raw_ms)
  raw_ms_info = raw_ms$ms_info %>%
    dplyr::filter(scan %in% raw_ms$scan_range)

  normalization_factors = zip_ms$peak_finder$peak_regions$normalization_factors
  final_scans = data.frame(scan = as.numeric(zip_ms$peak_finder$peak_regions$scan_level_arrays$Scan))

  freq_scans = names(zip_ms$peak_finder$peak_regions$frequency_point_regions$frequency)
  setdiff_raw_freq = length(setdiff(raw_ms_info$scan, freq_scans))
  setdiff_raw_norm = length(setdiff(raw_ms_info$scan, normalization_factors$scan))
  setdiff_freq_norm = length(setdiff(freq_scans, normalization_factors$scan))

  if (setdiff_raw_freq > 0) {
    mz_df_list = raw_ms$extract_raw_data()
    freq_coef2 = convert_scans2(mz_df_list, zip_ms$peak_finder$peak_regions$frequency_fit_description, zip_ms$peak_finder$peak_regions$mz_fit_description)
  } else {
    freq_coef2 = NULL
  }

  if (setdiff_freq_norm > 0) {
    out_diff = setdiff(freq_scans, normalization_factors$scan)
  } else {
    out_diff = NULL
  }

  zip_ms$cleanup()

  list(zip = in_zip,
       frequency = freq_coef2,
       normalization = out_diff)

}

checked_zip = purrr::map(all_zip, check_scans_removed)

saveRDS(checked_zip, file = "data/data_output/lung_data/lung_scandiffs_2022-05-11.rds")
