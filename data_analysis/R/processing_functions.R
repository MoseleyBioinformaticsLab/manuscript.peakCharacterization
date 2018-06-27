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
  char_ms$zip_ms$peak_finder$sample_id <- char_ms$zip_ms$id
  char_ms$zip_ms$cleanup()
  char_ms$zip_ms$peak_finder$add_sliding_regions()
  char_ms$zip_ms$peak_finder$add_tiled_regions()

  char_ms
}

