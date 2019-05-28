# common setup stuff
#
# read in the scans and generate the sliding and tiling windows
reading_scans_tile_windows <- function(in_file, pkg_desc){
  #use_file = file.path('data_analysis', 'data_input', in_file)
  use_file = in_file
  char_ms <- CharacterizeMSPeaks$new(use_file, peak_finder = PeakRegionFinder$new())
  char_ms$load_file()
  char_ms$filter_raw_scans()
  char_ms$zip_ms$peak_finder <- char_ms$peak_finder
  char_ms$zip_ms$peak_finder$add_data(char_ms$zip_ms$raw_ms)
  char_ms$zip_ms$raw_ms <- NULL # we drop the raw data to save us some memory
  char_ms$zip_ms$peak_finder$sample_id <- char_ms$zip_ms$id
  char_ms$zip_ms$cleanup() # don't want to leave stuff in the tmp directory
  char_ms$zip_ms$peak_finder$add_regions()

  list(char_obj = char_ms, pkg = pkg_desc)
}

reduce_removing_zero <- function(regions, point_regions, min_value = 0){
  regions <- FTMS.peakCharacterization:::count_overlaps(regions, point_regions)
  nz_counts <- regions@elementMetadata$nonzero_counts

  regions <- regions[nz_counts > min_value]
  IRanges::reduce(regions)
}

zero_normalization = function(peak_regions){
  intensity_measure = c("RawHeight", "Height")
  summary_function = median
  normalize_peaks = "both"
  scan_peaks <- peak_regions$scan_peaks

  all_scans = unique(unlist(purrr::map(scan_peaks, "scan")))
  normalization_factors <- data.frame(scan = all_scans, normalization = 0)

  normed_peaks <- internal_map$map_function(scan_peaks, normalize_scan_peaks, normalization_factors)

  normed_scan_cor <- purrr::map_dbl(normed_peaks, intensity_scan_correlation)
  normed_scan_cor[is.na(normed_scan_cor)] <- 0
  low_cor <- abs(normed_scan_cor) <= 0.5

  normed_raw <- normalize_raw_points(peak_regions$frequency_point_regions, normalization_factors)

  peak_regions$scan_peaks <- normed_peaks
  peak_regions$frequency_point_regions <- normed_raw
  peak_regions$is_normalized <- "both"
  peak_regions$normalization_factors <- normalization_factors

  normed_scan_cor <- data.frame(ScanCorrelation = normed_scan_cor,
                                HighCor = !low_cor)
  n_scans <- purrr::map_int(scan_peaks, calculate_number_of_scans)
  normed_scan_cor$HighScan <- n_scans >= quantile(n_scans, 0.9)
  normed_scan_cor$Ignore <- normed_scan_cor$HighCor & normed_scan_cor$HighScan
  peak_regions$scan_correlation <- normed_scan_cor
  peak_regions
}


group1_characterization = function(in_list){
  in_char = in_list$char_obj

  in_char$zip_ms$peak_finder$peak_regions$peak_regions =
    reduce_removing_zero(in_char$peak_finder$peak_regions$sliding_regions,
                         in_char$peak_finder$peak_regions$frequency_point_regions)
  in_char$zip_ms$peak_finder$peak_regions = zero_normalization(in_char$zip_ms$peak_finder$peak_regions)
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char
}

group2_characterization = function(in_list){
  in_char = in_list$char_obj

  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$peak_regions = zero_normalization(in_char$zip_ms$peak_finder$peak_regions)
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char
}

group3_characterization = function(in_list){
  in_char = in_list$char_obj
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()

  # do single pass normalization without considering intensity
  single_all_normalization = function(peak_regions){
    intensity_measure = c("RawHeight", "Height")
    summary_function = median
    normalize_peaks = "both"
    scan_peaks <- peak_regions$scan_peaks

    normalization_factors <- single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function,
                                                       min_ratio = 0)

    normed_peaks <- internal_map$map_function(scan_peaks, normalize_scan_peaks, normalization_factors)

    normed_scan_cor <- purrr::map_dbl(normed_peaks, intensity_scan_correlation)
    normed_scan_cor[is.na(normed_scan_cor)] <- 0
    low_cor <- abs(normed_scan_cor) <= 0.5

    normed_raw <- normalize_raw_points(peak_regions$frequency_point_regions, normalization_factors)

    peak_regions$scan_peaks <- normed_peaks
    peak_regions$frequency_point_regions <- normed_raw
    peak_regions$is_normalized <- "both"
    peak_regions$normalization_factors <- normalization_factors

    normed_scan_cor <- data.frame(ScanCorrelation = normed_scan_cor,
                                  HighCor = !low_cor)
    n_scans <- purrr::map_int(scan_peaks, calculate_number_of_scans)
    normed_scan_cor$HighScan <- n_scans >= quantile(n_scans, 0.9)
    normed_scan_cor$Ignore <- normed_scan_cor$HighCor & normed_scan_cor$HighScan
    peak_regions$scan_correlation <- normed_scan_cor
    peak_regions
  }

  in_char$zip_ms$peak_finder$peak_regions = single_all_normalization(in_char$zip_ms$peak_finder$peak_regions)
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char
}

group4_characterization = function(in_list){
  in_char = in_list$char_obj
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()

  # do single pass normalization with considering intensity
  single_high_normalization = function(peak_regions){
    intensity_measure = c("RawHeight", "Height")
    summary_function = median
    normalize_peaks = "both"
    scan_peaks <- peak_regions$scan_peaks

    normalization_factors <- single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function,
                                                       min_ratio = 0.7)

    normed_peaks <- internal_map$map_function(scan_peaks, normalize_scan_peaks, normalization_factors)

    normed_scan_cor <- purrr::map_dbl(normed_peaks, intensity_scan_correlation)
    normed_scan_cor[is.na(normed_scan_cor)] <- 0
    low_cor <- abs(normed_scan_cor) <= 0.5

    normed_raw <- normalize_raw_points(peak_regions$frequency_point_regions, normalization_factors)

    peak_regions$scan_peaks <- normed_peaks
    peak_regions$frequency_point_regions <- normed_raw
    peak_regions$is_normalized <- "both"
    peak_regions$normalization_factors <- normalization_factors

    normed_scan_cor <- data.frame(ScanCorrelation = normed_scan_cor,
                                  HighCor = !low_cor)
    n_scans <- purrr::map_int(scan_peaks, calculate_number_of_scans)
    normed_scan_cor$HighScan <- n_scans >= quantile(n_scans, 0.9)
    normed_scan_cor$Ignore <- normed_scan_cor$HighCor & normed_scan_cor$HighScan
    peak_regions$scan_correlation <- normed_scan_cor
    peak_regions
  }

  in_char$zip_ms$peak_finder$peak_regions = single_high_normalization(in_char$zip_ms$peak_finder$peak_regions)
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char
}

group5_characterization = function(in_list){
  in_char = in_list$char_obj
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()

  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char
}

group6_characterization = function(in_list){
  in_char = in_list$char_obj
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$remove_high_frequency_sd()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char
}

split_with_noise <- function(in_char){
  set_internal_map(furrr::future_map)
  out_char <- in_char$clone(deep = TRUE)
  out_char$peak_finder$peak_regions$peak_regions <-
    reduce_removing_zero(out_char$peak_finder$peak_regions$sliding_regions,
                         out_char$peak_finder$peak_regions$mz_point_regions)
  out_char$peak_finder$split_peak_regions()
  out_char$peak_finder$remove_double_peaks_in_scans()
  out_char
}

split_without_noise <- function(in_char){
  set_internal_map(furrr::future_map)
  out_char <- in_char$clone(deep = TRUE)
  out_char$peak_finder$reduce_sliding_regions()
  out_char$peak_finder$split_peak_regions()
  out_char$peak_finder$remove_double_peaks_in_scans()
  out_char
}

write_peaks_for_assignment <- function(in_regions, out_file){
  peak_list <- in_regions$summarize_peaks()
  cat(jsonlite::toJSON(peak_list, auto_unbox = TRUE, pretty = TRUE, digits = 8), file = out_file, sep = "\n")
}
