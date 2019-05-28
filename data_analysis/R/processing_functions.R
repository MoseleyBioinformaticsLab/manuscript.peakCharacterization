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

final_characterization = function(in_list){
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

reduce_removing_zero <- function(regions, point_regions, min_value = 0){
  regions <- FTMS.peakCharacterization:::count_overlaps(regions, point_regions)
  nz_counts <- regions@elementMetadata$nonzero_counts

  regions <- regions[nz_counts > min_value]
  IRanges::reduce(regions)
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

set_zero_normalization <- function(out_char){
  point_scans <- unique(out_char$peak_finder$peak_regions$mz_point_regions@elementMetadata$scan)
  norm_factors <- data.frame(scan = point_scans, normalization = 0)
  out_char$peak_finder$peak_regions$normalization_factors <- norm_factors
  out_char$peak_finder$peak_regions$is_normalized <- "both"
  out_char
}

run_characterization <- function(in_char, normalize = FALSE){
  set_internal_map(furrr::future_map)
  out_char <- in_char$clone(deep = TRUE)
  if (normalize) {
    out_char$peak_finder$normalize_data("both")
  } else {
    out_char <- set_zero_normalization(out_char)
  }
  out_char$peak_finder$find_peaks_in_regions()
  out_char$peak_finder$add_offsets()
  out_char$peak_finder$model_mzsd()
  out_char$peak_finder$model_heightsd()
}

write_peaks_for_assignment <- function(in_regions, out_file){
  peak_list <- in_regions$summarize_peaks()
  cat(jsonlite::toJSON(peak_list, auto_unbox = TRUE, pretty = TRUE, digits = 8), file = out_file, sep = "\n")
}
