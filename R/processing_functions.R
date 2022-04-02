# common setup stuff
#
# read in the scans and generate the sliding and tiling windows
reading_scans_tile_windows <- function(in_file, pkg){
  #use_file = file.path('data_analysis', 'data_input', in_file)
  use_file = in_file
  char_ms <- CharacterizeMSPeaks$new(use_file, peak_finder = PeakRegionFinder$new())
  char_ms$load_file()
  char_ms$peak_finder$start_time = Sys.time()
  char_ms$filter_raw_scans()
  char_ms$zip_ms$peak_finder <- char_ms$peak_finder
  char_ms$zip_ms$peak_finder$add_data(char_ms$zip_ms$raw_ms)
  char_ms$zip_ms$raw_ms <- NULL # we drop the raw data to save us some memory
  char_ms$zip_ms$peak_finder$sample_id <- char_ms$zip_ms$id
  char_ms$zip_ms$cleanup() # don't want to leave stuff in the tmp directory
  char_ms$zip_ms$peak_finder$add_regions()

  char_ms
}

reduce_removing_zero <- function(regions, point_regions_list, multiplier = 1.5, n_point_region = 10000){
  nz_counts = FTMS.peakCharacterization:::count_overlaps(regions, point_regions_list[[1]])
  n_region = seq(2, length(point_regions_list))
  for (iregion in n_region) {
    nz_counts_iregion = FTMS.peakCharacterization:::count_overlaps(regions, point_regions_list[[iregion]])
    nz_counts = nz_counts + nz_counts_iregion
  }

  regions = regions[nz_counts > 0]
  IRanges::reduce(regions)
}

find_signal_adj <- function(regions, point_regions_list, multiplier = 1.5, n_point_region = 2000, use_percentile = 0.98){
  nz_counts = FTMS.peakCharacterization:::count_overlaps(regions, point_regions_list[[1]])
  n_region = seq(2, length(point_regions_list))
  for (iregion in n_region) {
    nz_counts_iregion = FTMS.peakCharacterization:::count_overlaps(regions, point_regions_list[[iregion]])
    nz_counts = nz_counts + nz_counts_iregion
  }

  chunk_indices = seq(1, length(nz_counts), by = n_point_region)

  chunk_perc = purrr::map_dbl(chunk_indices, function(in_index){
    use_counts = nz_counts[seq(in_index, min(in_index + (n_point_region - 1), length(nz_counts)))]
    if (max(use_counts) > 0) {
      return(stats::quantile(use_counts, use_percentile))
    } else {
      return(0)
    }
  })

  cutoff_value = ceiling(median(chunk_perc) * multiplier)

  regions = regions[nz_counts > cutoff_value]
  IRanges::reduce(regions)
}


single_adjustable_normalization = function(peak_regions, min_ratio = 0){
  intensity_measure = c("RawHeight", "Height")
  summary_function = median
  normalize_peaks = "both"
  scan_peaks <- purrr::map(peak_regions$peak_region_list, "peaks")

  normalization_factors <- FTMS.peakCharacterization:::single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function,
                                                                                 min_ratio = min_ratio)
  named_norm = normalization_factors$normalization
  names(named_norm) = as.character(normalization_factors$scan)
  normed_list_regions = FTMS.peakCharacterization:::internal_map$map_function(peak_regions$peak_region_list, function(in_region){
    in_region$points = FTMS.peakCharacterization:::normalize_raw_points(in_region$points, named_norm)
    in_region$peaks = FTMS.peakCharacterization:::normalize_scan_peaks(in_region$peaks, normalization_factors)
    in_region
  })

  normed_raw <- FTMS.peakCharacterization:::normalize_raw_points(peak_regions$frequency_point_regions$frequency, named_norm)
  normed_peaks <- FTMS.peakCharacterization:::internal_map$map_function(scan_peaks, FTMS.peakCharacterization:::normalize_scan_peaks, normalization_factors)

  normed_scan_cor <- purrr::map_dbl(normed_peaks, FTMS.peakCharacterization:::intensity_scan_correlation)
  normed_scan_cor[is.na(normed_scan_cor)] <- 0
  low_cor <- abs(normed_scan_cor) <= 0.5


  peak_regions$peak_region_list = normed_list_regions
  peak_regions$frequency_point_regions$frequency <- normed_raw
  peak_regions$is_normalized <- "both"
  peak_regions$normalization_factors <- normalization_factors

  normed_scan_cor <- data.frame(ScanCorrelation = normed_scan_cor,
                                HighCor = !low_cor)
  n_scans <- purrr::map_int(scan_peaks, FTMS.peakCharacterization:::calculate_number_of_scans)
  peak_regions$scan_correlation <- normed_scan_cor
  peak_regions
}


# only remove zero points, no consideration of percentile, and no normalization (zero_normalization)
noperc_nonorm = function(use_char){
  in_char = use_char$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = TRUE
  in_char$zip_ms$peak_finder$progress = TRUE

  in_char$zip_ms$peak_finder$peak_regions$peak_regions =
    reduce_removing_zero(in_char$peak_finder$peak_regions$sliding_regions,
                          in_char$peak_finder$peak_regions$frequency_point_regions$frequency,
                         in_char$peak_finder$quantile_multiplier,
                         in_char$peak_finder$n_point_region)
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$indicate_high_frequency_sd()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char$zip_ms$peak_finder$model_mzsd()
  list(char_obj = in_char, processed = "noperc_nonorm")
}

# remove points based on the 99th percentile, no normalization (see zero_normalization)
perc99_nonorm = function(use_char){
  in_char = use_char$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = TRUE

  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$indicate_high_frequency_sd()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char$zip_ms$peak_finder$model_mzsd()
  list(char_obj = in_char, processed = "perc99_nonorm")
}

# remove points from 99th percentile, do single pass normalization without intensity cutoff
singlenorm = function(use_char){
  in_char = use_char$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = FALSE
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()

  # do single pass normalization without considering intensity

  in_char$zip_ms$peak_finder$peak_regions = single_adjustable_normalization(in_char$zip_ms$peak_finder$peak_regions, 0)
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$indicate_high_frequency_sd()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char$zip_ms$peak_finder$model_mzsd()
  list(char_obj = in_char, processed = "singlenorm")
}

# remove points below 99th percentile
# single pass normalization based on intensity
intsinglenorm = function(use_char){
  in_char = use_char$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = FALSE
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()

  # do single pass normalization with considering intensity
  in_char$zip_ms$peak_finder$peak_regions = single_adjustable_normalization(in_char$zip_ms$peak_finder$peak_regions, 0.7)
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$indicate_high_frequency_sd()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char$zip_ms$peak_finder$model_mzsd()
  list(char_obj = in_char, processed = "singlenorm_int")
}

# remove points based on 99th percentile
# proper two pass normalization
# no filtering on frequency sd
doublenorm = function(use_char){
  in_char = use_char$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = FALSE
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()

  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$indicate_high_frequency_sd()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char$zip_ms$peak_finder$model_mzsd()
  list(char_obj = in_char, processed = "doublenorm")
}

# remove points based on 99th percentile
# proper two pass normalization
# filter on frequency sd
filtersd = function(use_char){
  in_char = use_char$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = FALSE
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$indicate_high_frequency_sd()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char$zip_ms$peak_finder$model_mzsd()
  list(char_obj = in_char, processed = "filtersd")
}

filtersd98 = function(use_char){
  in_char = use_char$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = FALSE
  in_char$zip_ms$peak_finder$peak_regions$peak_regions =
    find_signal_98(in_char$peak_finder$peak_regions$sliding_regions,
                         in_char$peak_finder$peak_regions$frequency_point_regions$frequency,
                         in_char$peak_finder$quantile_multiplier,
                         in_char$peak_finder$n_point_region)
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$indicate_high_frequency_sd()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  in_char$zip_ms$peak_finder$model_mzsd()
  list(char_obj = in_char, processed = "filtersd98")
}

write_peaks_for_assignment <- function(in_data){
  in_char = in_data$char_obj
  processed_id = in_data$processed
  sample_id = in_char$zip_ms$peak_finder$sample_id
  out_id = dplyr::case_when(
  grepl("97Cpos", sample_id) ~ "97lipid",
  grepl("49Cpos", sample_id) ~ "49lipid",
  grepl("1_ECF", sample_id) ~ "1ecf",
  grepl("2_ECF", sample_id) ~ "2ecf")
  out_full = paste0(processed_id, "_", out_id, ".zip")
  write_loc = file.path("data/data_output", out_full)
  in_char$out_file = write_loc
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
  setwd("smirfe_related/smirfe_code/SMIRFE_assigner/")
  if (grepl("ecf", in_file)) {
    run_str = glue("python3.8 -m pipenv run python3 ./Main.py ../../smirfe_dbs/n15_1600_match_smirfemanuscript.db {zip_path} '_assigned.json' '[\"15N\"]' '[\"H\", \"Na\"]' '[1]'")
    system(run_str)
    setwd(curr_dir)
    out_value = return_file(assigned_file)
  } else {
    run_str = glue("python3.8 -m pipenv run python3 ./Main.py ../../smirfe_dbs/none_nosulfur_1600.db {zip_path} '_assigned.json' '[]' '[\"H\", \"Na\", \"K\", \"NH4\"]' '[1]'")
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
  evalue_cutoff = 0.1
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

  list(emfs = emfs$grouped_emf[possible_emfs$emf], summary_data = possible_emfs)

}

return_file = function(in_file){
  sha256 = digest::digest(algo = "sha256", file = in_file)
  list(file = in_file, sha256 = sha256)
}

normalization_factors = function(...){
  processed_data = list(...)

  get_factors = function(processed_obj){
    norm_df = processed_obj$char_obj$zip_ms$peak_finder$peak_regions$normalization_factors
    norm_df$processing = processed_obj$processed
    norm_df$sample_id = processed_obj$char_obj$zip_ms$peak_finder$sample_id
    norm_df
  }

  normalization_values = purrr::map_df(processed_data, get_factors)
  normalization_values

}

calc_rsd = function(in_char, processing){
  scan_peaks = 10^in_char$zip_ms$peak_finder$peak_regions$scan_level_arrays$Log10Height
  n_scan = ncol(scan_peaks)
  mean_peaks = rowMeans(scan_peaks, na.rm = TRUE)
  sd_peaks = apply(scan_peaks, 1, sd, na.rm = TRUE)
  n_peak = apply(scan_peaks, 1, function(.x){sum(!is.na(.x))})

  rsd_df = data.frame(mean = mean_peaks,
                      sd = sd_peaks,
                      rsd = sd_peaks / mean_peaks,
                      n = n_peak,
                      n_perc = n_peak / n_scan,
                      processed = processing,
                      sample = rename_samples( in_char$zip_ms$peak_finder$sample_id))
  rsd_df
}

single_rsd = function(char_list){
  rsd_df = calc_rsd(char_list$char_obj, char_list$processed)
  rsd_df
}

calc_rsd_msnbase = function(msnbase_list){
  scan_peaks = 10^msnbase_list$ScanLevel$Log10Height
  n_scan = ncol(scan_peaks)
  mean_peaks = rowMeans(scan_peaks, na.rm = TRUE)
  sd_peaks = apply(scan_peaks, 1, sd, na.rm = TRUE)
  n_peak = apply(scan_peaks, 1, function(.x){sum(!is.na(.x))})

  rsd_df = data.frame(mean = mean_peaks,
                      sd = sd_peaks,
                      rsd = sd_peaks / mean_peaks,
                      n = n_peak,
                      n_perc = n_peak / n_scan,
                      processed = dplyr::case_when(
                        grepl("pc", msnbase_list$Sample) ~ "msnbase_pc",
                        TRUE ~ "msnbase_only"
                      ),
                      sample = rename_samples(msnbase_list$Sample))
  rsd_df %>%
    dplyr::filter(n >= 3)
}

merge_rsd = function(...) {
  rsd_list = list(...)
  rsd_values = purrr::map_dfr(rsd_list, ~.x)
  rsd_values
}

summarize_rsd = function(rsd_df){
  get_mode = function(in_values){
    d_estimate = stats::density(in_values)
    y_max = which.max(d_estimate$y)
    d_estimate$x[y_max]
  }
  dplyr::group_by(rsd_df, processed, sample) %>%
    dplyr::summarize(mean = mean(rsd),
                     median = median(rsd),
                     mode = get_mode(rsd),
                     max = max(rsd))
}

split_regions = function(p99_nonorm_data){
  # loadd("method_perc99_nonorm_49lipid")
  # p99_nonorm_data = tar_read(method_perc99_nonorm_49lipid)


  pf1 = p99_nonorm_data$char_obj$zip_ms$peak_finder
  pf2 = p99_nonorm_data$char_obj$peak_finder

  pf3 = pf2$clone(deep = TRUE)
  pf3$reduce_sliding_regions()
  signal_regions = pf3$peak_regions$peak_regions

  frequency_point_regions = pf2$peak_regions$frequency_point_regions
  tiled_regions = pf1$peak_regions$tiled_regions
  min_scan = pf1$peak_regions$min_scan
  peak_method = pf1$peak_method
  min_points = pf1$min_points

  signal_list = as.list(split(signal_regions, seq(1, length(signal_regions))))
  min_scan2 = min(c(floor(min_scan / 2), 2))


  set.seed(1234)
  signal_list = signal_list[sample(length(signal_list), 500)]
  scan_regions_list = purrr::map(frequency_point_regions$frequency, function(in_points){

    points_list = purrr::map(signal_list, function(in_region){
      IRanges::subsetByOverlaps(in_points, in_region)
    })
    null_points = purrr::map_lgl(points_list, ~ length(.x) == 0)
    points_list[!null_points]
  })

  all_regions = vector("list", length(scan_regions_list))
  names(all_regions) = names(scan_regions_list)
  all_points = purrr::map(scan_regions_list, ~ names(.x)) %>% unlist(.) %>% unique(.)
  point_regions_list = vector("list", length(all_points))
  names(point_regions_list) = all_points

  pb = knitrProgressBar::progress_estimated(length(point_regions_list))

  for (ipoints in names(point_regions_list)) {
    tmp_scans = all_regions
    for (iscan in names(tmp_scans)) {
      tmp_scans[[iscan]] = scan_regions_list[[iscan]][[ipoints]]
    }
    nonzero_scans = purrr::map_dbl(tmp_scans, function(in_points){
      as.data.frame(in_points@elementMetadata) %>%
        dplyr::filter(intensity > 0) %>%
        dplyr::pull(scan) %>% unique(.) %>% length(.)
    })
    if (sum(nonzero_scans) >= min_scan2) {
      tiles = IRanges::subsetByOverlaps(tiled_regions, signal_list[[ipoints]])
      point_regions_list[[ipoints]] = list(
        points = tmp_scans[nonzero_scans > 0],
        tiles = tiles, region = signal_list[[ipoints]]
      )
    }
    knitrProgressBar::update_progress(pb)
  }


  not_null = purrr::map_lgl(point_regions_list, ~ !is.null(.x))
  point_regions_list = point_regions_list[not_null]
  point_regions_list = point_regions_list[sample(length(point_regions_list))]

  metadata = frequency_point_regions$metadata

  # split_data = purrr::imap(point_regions_list, function(.x, .y){
  #   message(.y)
  #   split_region_by_peaks(.x, peak_method = peak_method, min_points = min_points)
  # })
  split_data = FTMS.peakCharacterization:::internal_map$map_function(point_regions_list, FTMS.peakCharacterization:::split_region_by_peaks,
                                         peak_method = peak_method, min_points = min_points, metadata = metadata)

  null_regions = purrr::map_lgl(split_data, ~ is.null(.x[[1]]$points))
  split_data = split_data[!null_regions]

  n_regions = purrr::map_int(split_data, ~length(.x))

  possible_region = which(n_regions >= 2)
  n_each = purrr::imap_dfr(split_data[possible_region], function(.x, .y){
    purrr::map_dfr(seq(.x), function(in_x){
      data.frame(region = .y, sub_region = in_x,
                 n_peak = nrow(.x[[in_x]]$peaks))
    })
  })

  summary_region = n_each %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(min_scan = min(n_peak), n_region = n()) %>%
    dplyr::filter(min_scan > 40)

  use_list = point_regions_list[summary_region[["region"]]]

  use_list

}


rename_samples = function(sample_id){
  out_id = dplyr::case_when(
    grepl("97Cpos", sample_id) ~ "97lipid",
    grepl("49Cpos", sample_id) ~ "49lipid",
    grepl("1_ECF", sample_id) ~ "1ecf",
    grepl("2_ECF", sample_id) ~ "2ecf",
    grepl("1ecf", sample_id) ~ "1ecf",
    grepl("2ecf", sample_id) ~ "2ecf",
    grepl("97lipid", sample_id) ~ "97lipid",
    grepl("49lipid", sample_id) ~ "49lipid")
  out_id
}
