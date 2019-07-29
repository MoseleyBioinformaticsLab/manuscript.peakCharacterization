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

# only remove zero points, no consideration of percentile, and no normalization (zero_normalization)
group1_characterization = function(in_list){
  in_char = in_list$char_obj$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = TRUE

  in_char$zip_ms$peak_finder$peak_regions$peak_regions =
    reduce_removing_zero(in_char$peak_finder$peak_regions$sliding_regions,
                         in_char$peak_finder$peak_regions$frequency_point_regions)
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char
}

# remove points based on the 99th percentile, no normalization (see zero_normalization)
group2_characterization = function(in_list){
  in_char = in_list$char_obj$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = TRUE

  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char
}

# remove points from 99th percentile, do single pass normalization without intensity cutoff
group3_characterization = function(in_list){
  in_char = in_list$char_obj$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = FALSE
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()

  # do single pass normalization without considering intensity
  single_all_normalization = function(peak_regions){
    intensity_measure = c("RawHeight", "Height")
    summary_function = median
    normalize_peaks = "both"
    scan_peaks <- peak_regions$scan_peaks

    normalization_factors <- FTMS.peakCharacterization:::single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function,
                                                       min_ratio = 0)

    normed_peaks <- FTMS.peakCharacterization:::internal_map$map_function(scan_peaks, FTMS.peakCharacterization:::normalize_scan_peaks, normalization_factors)

    normed_scan_cor <- purrr::map_dbl(normed_peaks, FTMS.peakCharacterization:::intensity_scan_correlation)
    normed_scan_cor[is.na(normed_scan_cor)] <- 0
    low_cor <- abs(normed_scan_cor) <= 0.5

    normed_raw <- FTMS.peakCharacterization:::normalize_raw_points(peak_regions$frequency_point_regions, normalization_factors)

    peak_regions$scan_peaks <- normed_peaks
    peak_regions$frequency_point_regions <- normed_raw
    peak_regions$is_normalized <- "both"
    peak_regions$normalization_factors <- normalization_factors

    normed_scan_cor <- data.frame(ScanCorrelation = normed_scan_cor,
                                  HighCor = !low_cor)
    n_scans <- purrr::map_int(scan_peaks, FTMS.peakCharacterization:::calculate_number_of_scans)
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

# remove points below 99th percentile
# single pass normalization based on intensity
group4_characterization = function(in_list){
  in_char = in_list$char_obj$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = FALSE
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()

  # do single pass normalization with considering intensity
  single_high_normalization = function(peak_regions){
    intensity_measure = c("RawHeight", "Height")
    summary_function = median
    normalize_peaks = "both"
    scan_peaks <- peak_regions$scan_peaks

    normalization_factors <- FTMS.peakCharacterization:::single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function,
                                                       min_ratio = 0.7)

    normed_peaks <- FTMS.peakCharacterization:::internal_map$map_function(scan_peaks, FTMS.peakCharacterization:::normalize_scan_peaks, normalization_factors)

    normed_scan_cor <- purrr::map_dbl(normed_peaks, FTMS.peakCharacterization:::intensity_scan_correlation)
    normed_scan_cor[is.na(normed_scan_cor)] <- 0
    low_cor <- abs(normed_scan_cor) <= 0.5

    normed_raw <- FTMS.peakCharacterization:::normalize_raw_points(peak_regions$frequency_point_regions, normalization_factors)

    peak_regions$scan_peaks <- normed_peaks
    peak_regions$frequency_point_regions <- normed_raw
    peak_regions$is_normalized <- "both"
    peak_regions$normalization_factors <- normalization_factors

    normed_scan_cor <- data.frame(ScanCorrelation = normed_scan_cor,
                                  HighCor = !low_cor)
    n_scans <- purrr::map_int(scan_peaks, FTMS.peakCharacterization:::calculate_number_of_scans)
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

# remove points based on 99th percentile
# proper two pass normalization
# no filtering on frequency sd
group5_characterization = function(in_list){
  in_char = in_list$char_obj$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = FALSE
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()

  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char
}

# remove points based on 99th percentile
# proper two pass normalization
# filter on frequency sd
group6_characterization = function(in_list){
  in_char = in_list$char_obj$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = FALSE
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

write_peaks_for_assignment <- function(in_char){
  in_char$out_file = file.path(write_loc, paste0(in_char$zip_ms$peak_finder$sample_id, ".zip"))
  in_char$zip_ms$out_file = in_char$out_file
  dir.create(in_char$zip_ms$temp_directory, recursive = TRUE)
  in_char$zip_ms$peak_finder$start_time = Sys.time()
  in_char$zip_ms$peak_finder$stop_time = Sys.time()
  in_char$summarize()
  in_char$save_peaks()
  in_char$write_zip()
  in_char$zip_ms$cleanup()
  return_file(in_char$out_file)
}

assign_files = function(in_zip){
  in_file = in_zip$file
  zip_path = normalizePath(in_file)
  curr_dir = getwd()
  assigned_file = paste0(in_file, "_assigned.json")

  if (grepl("ECF", in_file)) {
    setwd("~/Projects/work/SMIRFE/SMIRFE_assigner/")
    run_str = glue("python3 ./Main.py /mlab/scratch/cesb_data/smirfe_dbs/n15_1600.db {zip_path} '_assigned.json' '[\"15N\"]' '[\"H\", \"Na\", \"K\", \"NH4\"]' '[1]'")
    system(run_str)
    setwd(curr_dir)
    out_value = return_file(assigned_file)
  } else {
    setwd("~/Projects/work/SMIRFE/SMIRFE_assigner/")
    run_str = glue("python3 ./Main.py /mlab/scratch/cesb_data/smirfe_dbs/none_1600.db {zip_path} '_assigned.json' '[]' '[\"H\", \"Na\", \"K\", \"NH4\"]' '[1]'")
    system(run_str)
    setwd(curr_dir)
    out_value = return_file(assigned_file)
  }

  out_value
}

read_assignments = function(in_assign){
  if (!is.null(in_assign)) {
    sample_assignments = read_smirfe_assignment(in_assign$file)
  } else {
    sample_assignments = NULL
  }
  sample_assignments
}

find_interesting_peaks = function(assigned_data){
  evalue_cutoff = 0.5
  remove_elements = "S"
  tmp_assign = dplyr::filter(assigned_data$assignments, !grepl(remove_elements, complete_EMF))
  emfs = get_sample_emfs(tmp_assign, assigned_data$sample, evalue_cutoff = evalue_cutoff)

  peak_nscan = purrr::map_int(assigned_data$scan_level$ObservedFrequency, ~ sum(!is.na(.x)))

  organized_data = purrr::imap_dfr(emfs$grouped_emf, function(in_emf, emf_id){
    tmp_out = data.frame(emf = emf_id,
               clique_size = in_emf$clique_size,
               n_peak = length(in_emf$Sample_Peak),
               e_value = in_emf$min_e_value,
               stringsAsFactors = FALSE)
    tmp_out$peaks = list(in_emf$Sample_Peak)
    tmp_out$scans = list(peak_nscan[in_emf$Sample_Peak])
    tmp_out$min_scan = min(peak_nscan[in_emf$Sample_Peak])
    tmp_out
  })
  max_scan = max(unlist(purrr::map(organized_data$scans, ~ .x)))
  possible_emfs = dplyr::filter(organized_data, clique_size > 3, min_scan < (0.5 * max_scan))
  possible_emfs

}

return_file = function(in_file){
  sha256 = system2("sha256sum", args = in_file, stdout = TRUE)
  list(file = in_file, sha256 = sha256)
}
