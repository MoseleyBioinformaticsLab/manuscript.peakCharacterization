# common setup stuff
#
# average the scans and generate the sliding and tiling windows
average_scans_tile_windows <- function(in_file, pkg_description){
  char_ms <- CharacterizeMSPeaks$new(in_file, peak_finder = PeakRegionFinder$new())
  char_ms$load_file()
  set_internal_map(furrr::future_map)
  char_ms$filter_raw_scans()
  char_ms$zip_ms$peak_finder <- char_ms$peak_finder
  char_ms$zip_ms$peak_finder$add_data(char_ms$zip_ms$raw_ms)
  char_ms$zip_ms$raw_ms <- NULL # we drop the raw data to save us some memory
  char_ms$zip_ms$peak_finder$sample_id <- char_ms$zip_ms$id
  char_ms$zip_ms$cleanup() # don't want to leave stuff in the tmp directory
  char_ms$zip_ms$peak_finder$add_sliding_regions()
  char_ms$zip_ms$peak_finder$add_tiled_regions()

  char_ms
}

reduce_removing_zero <- function(regions, point_regions, min_value = 0){
  regions <- FTMS.peakCharacterization:::count_overlaps(regions, point_regions)
  nz_counts <- regions@elementMetadata$nonzero_counts

  regions <- regions[nz_counts > min_value]
  IRanges::reduce(regions)
}

split_with_noise <- function(in_char){
  out_char <- in_char$clone(deep = TRUE)
  out_char$peak_finder$peak_regions$peak_regions <-
    reduce_removing_zero(out_char$peak_finder$peak_regions$sliding_regions,
                         out_char$peak_finder$peak_regions$mz_point_regions)
  out_char$peak_finder$split_peak_regions()
  out_char
}

split_without_noise <- function(in_char){
  out_char <- in_char$clone(deep = TRUE)
  out_char$peak_finder$reduce_sliding_regions()
  out_char$peak_dinfer$split_peak_regions()
  out_char
}

run_characterization <- function(in_char, normalize = FALSE){
  out_char <- in_char$clone(deep = TRUE)
  if (normalize) {
    out_char$peak_finder$normalize_data()
  }
  out_char$peak_finder$find_peaks_in_regions()
  out_char$peak_finder$add_offsets()
  out_char$peak_finder$model_mzsd()
  out_char$peak_finder$model_heightsd()
}


