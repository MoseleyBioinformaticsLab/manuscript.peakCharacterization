#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param coefficients_data
#' @return
#' @author rmflight
#' @export
coefficients_analysis <- function(coefficients_data) {

  n_freq = purrr::map_dbl(coefficients_data, function(.x){
    length(.x$frequency$bad)
  })

  example_freq = coefficients_data[[(which(n_freq == 5)[1])]]

  df_freq = data.frame(coefficient = example_freq$frequency$coefficients,
                       sample = gsub(".zip$", "", basename(example_freq$zip)))
  df_freq = df_freq %>%
    dplyr::mutate(outlier = coefficient %in% example_freq$frequency$bad)

  n_norm = purrr::map_dbl(coefficients_data, function(.x){
    length(.x$normalization)
  })

  list(n_coef_out = n_freq,
       example = df_freq)

}

resolution_analysis = function(zip_sample){
  # zip_sample = "data/data_output/149Cpos.zip"
  zip_ms = ZipMS$new(zip_sample)
  zip_ms$load_raw()
  zip_ms$load_peak_finder()

  raw_info = zip_ms$raw_ms$ms_info
  raw_info$scan = as.character(raw_info$scan)

  zip_ms$raw_ms$set_scans()

  mz_df_list = zip_ms$raw_ms$extract_raw_data()
  frequency_fit_description = zip_ms$peak_finder$peak_regions$frequency_fit_description
  mz_fit_description = zip_ms$peak_finder$peak_regions$mz_fit_description
  ... = NULL

  if (is.null(names(mz_df_list))) {
    names(mz_df_list) = purrr::map_chr(mz_df_list, ~ .x$scan[1])
  }

  mz_frequency = FTMS.peakCharacterization:::internal_map$map_function(mz_df_list, function(in_scan){
    #message(scan_name)
    out_scan = FTMS.peakCharacterization:::convert_mz_frequency(in_scan, ...)
    out_scan
  })
  #str(mz_frequency[[1]])
  #log_message("converted frequencies")
  frequency_fits = FTMS.peakCharacterization:::internal_map$map_function(mz_frequency, function(in_freq){
    use_peaks = in_freq$convertable
    tmp_fit = FTMS.peakCharacterization:::fit_exponentials(in_freq$mean_mz[use_peaks], in_freq$mean_frequency[use_peaks], frequency_fit_description)
    tmp_fit$scan = in_freq[1, "scan"]
    tmp_fit
  })
  #str(frequency_fits[[1]])
  #log_message("fit frequencies")

  mz_fits = FTMS.peakCharacterization:::internal_map$map_function(mz_frequency, function(in_freq){
    use_peaks = in_freq$convertable
    tmp_fit = FTMS.peakCharacterization:::fit_exponentials(in_freq$mean_frequency[use_peaks], in_freq$mean_mz[use_peaks], mz_fit_description)
    tmp_fit$scan = in_freq[1, "scan"]
    tmp_fit
  })
  #str(mz_fits[[1]])
  #log_message("fit mz")

  frequency_coefficients = purrr::map_df(frequency_fits, function(.x){
    tmp_df = as.data.frame(matrix(.x$coefficients, nrow = 1))
    tmp_df$scan = .x$scan
    tmp_df
  })
  #str(frequency_coefficients)
  #log_message("extracted coefficients")

  mz_coefficients = purrr::map_df(mz_fits, function(.x){
    tmp_df = as.data.frame(matrix(.x$coefficients, nrow = 1))
    tmp_df$scan = .x$scan
    tmp_df
  })
  #str(mz_coefficients)

  first_slope = min(which(frequency_fit_description != 0))
  bad_coefficients = boxplot.stats(frequency_coefficients[[first_slope]])$out

  frequency_coefficients = frequency_coefficients[!(frequency_coefficients[[first_slope]] %in% bad_coefficients), ]
  frequency_coefficients$scan = as.character(frequency_coefficients$scan)
  frequency_coefficients = frequency_coefficients %>%
    dplyr::mutate(sqrt = dplyr::case_when(
      dplyr::between(V2, 29800000, 29900000) ~ "2.9e7",
      dplyr::between(V2, 7450000, 7460000) ~ "7.4e6"
    ))
  raw_info = dplyr::left_join(raw_info, frequency_coefficients, by = "scan")
  zip_ms$cleanup()
  raw_info
}
