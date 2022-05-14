model_rtime_filter = function(scan_times, min_time_difference = 4, rtime_limit = 7.5*60){

  low_rtime = scan_times$rtime <= rtime_limit

  high_lag <- scan_times$rtime_lag >= min_time_difference
  high_lag[is.na(high_lag)] <- TRUE
  high_lead <- scan_times$rtime_lead >= min_time_difference
  high_lead[is.na(high_lead)] <- TRUE

  na_lead_high_lag <- is.na(scan_times$rtime_lead) & high_lag
  na_lag_high_lead <- is.na(scan_times$rtime_lag) & high_lead

  keep_scans <- low_rtime & (na_lead_high_lag | high_lag) & (na_lag_high_lead | high_lead)
  keep_scans
}

model_time_filter = function(scan_times, min_time_difference = 4){

  high_lag <- scan_times$rtime_lag >= min_time_difference
  high_lag[is.na(high_lag)] <- TRUE
  high_lead <- scan_times$rtime_lead >= min_time_difference
  high_lead[is.na(high_lead)] <- TRUE

  na_lead_high_lag <- is.na(scan_times$rtime_lead) & high_lag
  na_lag_high_lead <- is.na(scan_times$rtime_lag) & high_lead

  keep_scans <- (na_lead_high_lag | high_lag) & (na_lag_high_lead | high_lead)
  keep_scans
}


run_frequency_models = function(mzml_file){
  mod_1 = c("a.freq" = 0, "y.freq" = -1/2)
  mod_2 = c("a.freq" = 0, "x.freq" = -1, "y.freq" = -1/2)
  mod_3 = c("a.freq" = 0, "x.freq" = -1, "y.freq" = -1/2, "z.freq" = -1/3)
  mod_4 = c("a.freq" = 0, "y.freq" = -1/2, "z.freq" = -1/3)

  list_models = list(mod1 = mod_1,
                     mod2 = mod_2,
                     mod3 = mod_3,
                     mod4 = mod_4)

  if (grepl("pos|lipid", mzml_file)) {
    use_filter = model_rtime_filter
  } else {
    use_filter = model_time_filter
  }


  sc_raw = ScanCentricPeakCharacterization::SCRaw$new(mzml_file)
  sc_raw$extract_raw_data()

  sc_raw$frequency_fit_description = mod_3
  sc_raw$predict_frequency()
  mod_3_freq = sc_raw$scan_info
  use_scans = use_filter(mod_3_freq)
  sc_raw$raw_df_data = sc_raw$raw_df_data[use_scans]
  sc_raw$scan_info = sc_raw$scan_info[use_scans, ]
  for (iname in names(mod_3)) {
    sc_raw$scan_info[[iname]] = NULL
  }

  out_models = purrr::map(list_models, function(.x){
    sc_raw$frequency_fit_description = .x
    sc_raw$predict_frequency()
    tmp_out = list(scans = sc_raw$raw_df_data,
                   freq_data = sc_raw$scan_info)
    for (iname in names(.x)) {
      sc_raw$scan_info[[iname]] = NULL
    }
    tmp_out
  })

  all_mads = purrr::imap_dfr(out_models, function(.x, .y){
    tmp_info = .x$freq_data %>%
      dplyr::select(scan, median, mad) %>%
      dplyr::mutate(model = .y)
    tmp_info
  })

  all_coef = names(list_models$mod3)
  coef_df = as.data.frame(matrix(NA, nrow = 1, ncol = length(all_coef)))
  names(coef_df) = all_coef

  get_coef = function(curr_coef, coef_df){
    match_coef = names(coef_df)[names(coef_df) %in% names(curr_coef)]
    tmp_coef = coef_df
    for (icoef in match_coef) {
      tmp_coef[[icoef]] = curr_coef[[icoef]][1]
    }
    tmp_coef
  }

  model_coefficients = purrr::imap_dfr(out_models, function(in_model, id){
    out_coef = get_coef(in_model$freq_data, coef_df)
    out_coef$model = id
    out_coef
  })

  sample_id = rename_samples(basename(mzml_file))
  list(sample = sample_id,
       mad_data = all_mads,
       model_data = model_coefficients,
       models = list_models)
}

combine_frequency_models = function(...){
  sample_list = list(...)
  sample_list
}
