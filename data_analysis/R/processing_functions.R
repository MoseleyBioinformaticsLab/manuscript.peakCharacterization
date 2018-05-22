raw_peakpicking <- function(zip_file, pkg_description){
  zip_data <- zip_ms(zip_file)
  zip_data$load_raw()
  zip_data$cleanup()
  zip_data$peak_finder <- PeakFinder$new()

  zip_data$peak_finder$raw_data <- zip_data$raw_ms

  zip_data$peak_finder$create_multi_scan_peaklist()
  zip_data
}

# identity offset correction, i.e. do no offset correction
# because we want to keep lots of other information, we are just replacing
# the correction factor with 0, and subtracting that out
identity_correct_offset_function <- function(master_peak_list, multi_scan_peaklist, min_scan = 0.1){

  n_col <- ncol(master_peak_list$scan_mz)
  n_min_scan <- SIRM.FTMS.peakCharacterization:::correctly_round_numbers(n_col, min_scan)
  correspond_peaks <- master_peak_list$count_notna() >= n_min_scan

  scan_indices <- master_peak_list$scan_indices

  multi_scan_peaklist_corrected <- multi_scan_peaklist$clone(deep = TRUE)

  offset_models <- vector(mode = "list", length = length(scan_indices))

  for (iscan in seq(1, n_col)) {
    scan_offset <- master_peak_list$scan_mz[correspond_peaks, iscan] -
      master_peak_list$master[correspond_peaks]
    offset_model <- master_peak_list$offset_fit_function(master_peak_list$master[correspond_peaks], scan_offset)

    use_scan <- scan_indices[iscan]
    tmp_peaks <- multi_scan_peaklist_corrected$peak_list_by_scans[[use_scan]]$peak_list
    use_mz <- tmp_peaks$ObservedMZ
    offset_predict <- 0
    corrected_mz <- use_mz - offset_predict
    tmp_peaks$ObservedMZ <- corrected_mz
    multi_scan_peaklist_corrected$peak_list_by_scans[[use_scan]]$peak_list <- tmp_peaks
    offset_model$scan_index <- scan_indices[iscan]
    offset_models[[iscan]] <- offset_model
  }

  list(multi_scan_peaklist = multi_scan_peaklist_corrected,
       models = offset_models)
}

keep_noise_no_offset_correspondence <- function(zip_data){
  zip_data2 <- zip_data$clone(deep = TRUE)
  zip_data2$raw_ms <- NULL
  zip_data2$peak_finder$multi_scan_peaklist$peak_list_by_scans <- purrr::map(zip_data2$peak_finder$multi_scan_peaklist$peak_list_by_scans,
                                                                             function(x){
                                                                               x$peak_list$not_noise <- TRUE
                                                                               x
                                                                             })
  zip_data2$peak_finder$keep_all_master_peak_lists <- TRUE
  zip_data2$peak_finder$keep_intermediates <- TRUE
  zip_data$peak_finder$offset_correction_function <- identity_correct_offset_function
  zip_data2$peak_finder$create_correspondent_peaks()
  zip_data2
}

remove_noise_no_offset_correspondence <- function(zip_data){
  zip_data2 <- zip_data$clone(deep = TRUE)
  zip_data2$raw_ms <- NULL
  zip_data$peak_finder$offset_correction_function <- identity_correct_offset_function
  zip_data2$peak_finder$keep_all_master_peak_lists <- TRUE
  zip_data2$peak_finder$keep_intermediates <- TRUE
  zip_data2$peak_finder$create_correspondent_peaks()
  zip_data2
}

remove_noise_do_offset_correspondence <- function(zip_data){
  zip_data2 <- zip_data$clone(deep = TRUE)
  zip_data2$raw_ms <- NULL
  zip_data2$peak_finder$keep_all_master_peak_lists <- TRUE
  zip_data2$peak_finder$keep_intermediates <- TRUE
  zip_data2$peak_finder$create_correspondent_peaks()
  zip_data2
}
