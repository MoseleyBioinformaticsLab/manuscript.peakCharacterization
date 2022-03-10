msnbase_centroid = function(mzml_file, char_list = NULL){


  if (!is.null(char_list)) {
    # get the scans from the normalization object,
    # and then subset the profiled scans using it
  }

  mzml_prof = MSnbase::readMSData(mzml_file, msLevel. = 1, centroided. = FALSE)
  mzml_info = get_ms_info(mzml_prof)
  mzml_info_hitime = msnbase_time_filter(mzml_info, 4)

  mzml_prof2 = mzml_prof[mzml_info_hitime$scan]

  all_scans_cent = mzml_prof2 %>%
    MSnbase::pickPeaks()
  all_scans_cent_mz = mz(all_scans_cent)
  all_scans_cent_intensity = intensity(all_scans_cent)

  all_scans_data = purrr::map_df(seq_len(length(all_scans_cent_mz)), function(in_scan){
    data.frame(mz = all_scans_cent_mz[[in_scan]],
               intensity = all_scans_cent_intensity[[in_scan]],
               scan = in_scan)
  })

  comb_prof = MSnbase::combineSpectra(mzml_prof2, method = meanMzInts, mzd = 0, ppm = 1)

  comb_cent = comb_prof %>%
    MSnbase::pickPeaks()
  comb_cent_data = data.frame(mz = mz(comb_cent)[[1]],
                              intensity = intensity(comb_cent)[[1]])

  list(scanlevel = all_scans_data, comb = comb_cent_data)
}

msnbase_time_filter <- function(scan_times, min_time_difference = 4){
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
