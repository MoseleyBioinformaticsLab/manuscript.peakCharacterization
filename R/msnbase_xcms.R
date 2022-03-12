msnbase_centroid = function(mzml_file, char_list = NULL){
  mzml_base = basename(mzml_file)
  sample_id = dplyr::case_when(
    grepl("97Cpos", mzml_base) ~ "97lipid",
    grepl("49Cpos", mzml_base) ~ "49lipid",
    grepl("1_ECF", mzml_base) ~ "1ecf",
    grepl("2_ECF", mzml_base) ~ "2ecf"
  )

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

  list(scanlevel = all_scans_data, comb = comb_cent_data,
       sample_id = sample_id)
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

msnbase_match_scan_combined = function(msnbase_data){
  comb_data = msnbase_data$comb %>%
    dplyr::mutate(ppm1 = mz * 2e-6) %>%
    dplyr::mutate(low = mz - ppm1,
                  hi = mz + ppm1) %>%
  dplyr::mutate(PeakID = as.character(seq(1, nrow(.))))

  scan_data = msnbase_data$scanlevel
  n_scan = length(unique(scan_data$scan))
  scan_mz = scan_intensity = matrix(NA, nrow = nrow(comb_data), ncol = n_scan)
  rownames(scan_mz) = rownames(scan_intensity) = comb_data$PeakID
  colnames(scan_mz) = colnames(scan_intensity) = unique(scan_data$scan)

  for (ipeak in seq_len(nrow(comb_data))) {
    scan_between = dplyr::between(scan_data$mz, comb_data[ipeak, "low"], comb_data[ipeak, "hi"])
    tmp_scan = scan_data[scan_between, ]
    scan_mz[comb_data[ipeak, "PeakID"], tmp_scan$scan] = tmp_scan$mz
    scan_intensity[comb_data[ipeak, "PeakID"], tmp_scan$scan] = tmp_scan$intensity
  }

  n_scan_peak = apply(scan_mz, 1, function(.x){
    sum(!is.na(.x))
  })
  comb_data = comb_data[n_scan_peak > 0, ]
  scan_mz = scan_mz[n_scan_peak > 0, ]
  scan_intensity = scan_intensity[n_scan_peak > 0, ]
  n_scan_peak = n_scan_peak[n_scan_peak > 0]
  new_id = paste0("msnbase_", msnbase_data$sample_id)
  new_data = comb_data %>%
    dplyr::transmute(PeakID = PeakID,
                     Height = intensity,
                     ObservedMZ = mz,
                     Log10Height = log10(intensity),
                     CorrectedLog10Height = log10(intensity),
                     NScan = n_scan_peak)
  scan_level_list = list(Log10Height = log10(scan_intensity),
                         CorrectedLog10Height = log10(scan_intensity),
                         ObservedMZ = scan_mz,
                         Scan = colnames(scan_mz),
                         PeakID = rownames(scan_mz)
                         )
  list(TIC = sum(new_data$Height),
       Sample = new_id,
       Peaks = new_data,
       ScanLevel = scan_level_list)
}

msnbase_match_pc = function(msnbase_data, char_list){
  peak_data = char_list$char_obj$zip_ms$peak_finder$peak_regions$peak_data

  msnbase_peaks = msnbase_data$comb
  msnbase_peaks$PeakID = as.character(seq(1, nrow(msnbase_peaks)))
  msnbase_peaks$PeakCharID = "NA"
  msnbase_scans = msnbase_data$scanlevel

  n_scan = length(unique(msnbase_scans$scan))

  model_matrix = matrix(NA, nrow = 1, ncol = n_scan)
  colnames(model_matrix) = unique(msnbase_scans$scan)

  matched_peak_list = purrr::map(seq_len(nrow(peak_data)), function(in_row){
    found_peaks = dplyr::between(msnbase_peaks$mz,
                                 peak_data[in_row, "Start"],
                                 peak_data[in_row, "Stop"])
    found_scan = dplyr::between(msnbase_scans$mz,
                                peak_data[in_row, "Start"],
                                peak_data[in_row, "Stop"])

    # if we find more than a single,
    #  OR
    # find nothing, don't return anything
    if ((sum(found_peaks) > 1) || (sum(found_peaks) == 0) || (sum(found_scan) == 0)) {
      return(NULL)
    }
    # check for multiple peaks in the same scans
    scan_df = msnbase_scans[found_scan, , drop = FALSE]
    dup_scan = scan_df$scan[duplicated(scan_df$scan)]
    if (length(dup_scan) > 0) {
      scan_df = scan_df %>%
        dplyr::filter(!(scan %in% dup_scan))
    }
    if (nrow(scan_df) == 0) {
      return(NULL)
    }

    tmp_mz = tmp_intensity = model_matrix
    rownames(tmp_mz) = rownames(tmp_intensity) = msnbase_peaks$PeakID[found_peaks]
    tmp_mz[, scan_df$scan] = scan_df$mz
    tmp_intensity[, scan_df$scan] = scan_df$intensity
    tmp_msnbase = msnbase_peaks[found_peaks, , drop = FALSE]
    tmp_msnbase$PeakCharID = peak_data[in_row, "PeakID"]
    list(peaks = tmp_msnbase,
         mz = tmp_mz,
         intensity = tmp_intensity)
  })

  has_peaks = purrr::map_lgl(matched_peak_list, ~ !is.null(.x))
  matched_peak_list = matched_peak_list[has_peaks]

  scan_mz = do.call(rbind, purrr::map(matched_peak_list, ~ .x$mz))
  scan_intensity = do.call(rbind, purrr::map(matched_peak_list, ~ .x$intensity))
  n_scan_peak = apply(scan_mz, 1, function(.x){
    sum(!is.na(.x))
  })

  cent_data = purrr::map_df(matched_peak_list, ~ .x$peaks)

  new_data = cent_data %>%
    dplyr::transmute(PeakID = PeakID,
                     Height = intensity,
                     ObservedMZ = mz,
                     Log10Height = log10(intensity),
                     CorrectedLog10Height = log10(intensity),
                     NScan = n_scan_peak,
                     PeakCharID = PeakCharID)
  scan_level_list = list(Log10Height = log10(scan_intensity),
                         CorrectedLog10Height = log10(scan_intensity),
                         ObservedMZ = scan_mz,
                         Scan = colnames(scan_mz),
                         PeakID = rownames(scan_mz)
  )
  new_id = paste0("pc_msnbase_", msnbase_data$sample_id)
  list(TIC = sum(new_data$Height),
       Sample = new_id,
       Peaks = new_data,
       ScanLevel = scan_level_list)
}

msnbase_zip = function(msnbase_list){
  out_list = list(peak_list.json = msnbase_list)
  zip_file = paste0(msnbase_list$Sample, ".zip")
  zip_loc = file.path("data/data_output", zip_file)
  zip_save = FTMS.peakCharacterization::lists_2_json(out_list, zip_file = zip_loc)
  return_file(zip_save)
}
