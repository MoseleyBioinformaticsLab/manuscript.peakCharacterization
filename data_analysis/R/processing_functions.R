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
noperc_nonorm = function(in_list){
  in_char = in_list$char_obj$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = TRUE
  in_char$zip_ms$peak_finder$progress = TRUE

  in_char$zip_ms$peak_finder$peak_regions$peak_regions =
    reduce_removing_zero(in_char$peak_finder$peak_regions$sliding_regions,
                         in_char$peak_finder$peak_regions$frequency_point_regions)
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  list(char_obj = in_char, processed = "noperc_nonorm")
}

# remove points based on the 99th percentile, no normalization (see zero_normalization)
perc99_nonorm = function(in_list){
  in_char = in_list$char_obj$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = TRUE

  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()
  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  list(char_obj = in_char, processed = "perc99_nonorm")
}

# remove points from 99th percentile, do single pass normalization without intensity cutoff
singlenorm = function(in_list){
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
    scan_peaks <- purrr::map(peak_regions$peak_region_list, "peaks")

    normalization_factors <- FTMS.peakCharacterization:::single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function,
                                                       min_ratio = 0)

    normed_list_regions = FTMS.peakCharacterization:::internal_map$map_function(peak_regions$peak_region_list, function(in_region){
      in_region$points = FTMS.peakCharacterization:::normalize_raw_points(in_region$points, normalization_factors)
      in_region$peaks = FTMS.peakCharacterization:::normalize_scan_peaks(in_region$peaks, normalization_factors)
      in_region
    })

    normed_peaks <- FTMS.peakCharacterization:::internal_map$map_function(scan_peaks, FTMS.peakCharacterization:::normalize_scan_peaks, normalization_factors)

    normed_scan_cor <- purrr::map_dbl(normed_peaks, FTMS.peakCharacterization:::intensity_scan_correlation)
    normed_scan_cor[is.na(normed_scan_cor)] <- 0
    low_cor <- abs(normed_scan_cor) <= 0.5

    normed_raw <- FTMS.peakCharacterization:::normalize_raw_points(peak_regions$frequency_point_regions, normalization_factors)

    peak_regions$peak_region_list = normed_list_regions
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
  list(char_obj = in_char, processed = "singlenorm")
}

# remove points below 99th percentile
# single pass normalization based on intensity
singlenorm_int = function(in_list){
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
    scan_peaks <- purrr::map(peak_regions$peak_region_list, "peaks")

    normalization_factors <- FTMS.peakCharacterization:::single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function,
                                                       min_ratio = 0.7)

    normed_list_regions = FTMS.peakCharacterization:::internal_map$map_function(peak_regions$peak_region_list, function(in_region){
      in_region$points = FTMS.peakCharacterization:::normalize_raw_points(in_region$points, normalization_factors)
      in_region$peaks = FTMS.peakCharacterization:::normalize_scan_peaks(in_region$peaks, normalization_factors)
      in_region
    })

    normed_peaks <- FTMS.peakCharacterization:::internal_map$map_function(scan_peaks, FTMS.peakCharacterization:::normalize_scan_peaks, normalization_factors)

    normed_scan_cor <- purrr::map_dbl(normed_peaks, FTMS.peakCharacterization:::intensity_scan_correlation)
    normed_scan_cor[is.na(normed_scan_cor)] <- 0
    low_cor <- abs(normed_scan_cor) <= 0.5

    normed_raw <- FTMS.peakCharacterization:::normalize_raw_points(peak_regions$frequency_point_regions, normalization_factors)

    peak_regions$peak_region_list = normed_list_regions
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
  list(char_obj = in_char, processed = "singlenorm_int")
}

# remove points based on 99th percentile
# proper two pass normalization
# no filtering on frequency sd
doublenorm = function(in_list){
  in_char = in_list$char_obj$clone(deep = TRUE)
  in_char$zip_ms$peak_finder$zero_normalization = FALSE
  in_char$zip_ms$peak_finder$reduce_sliding_regions()
  in_char$zip_ms$peak_finder$split_peak_regions()
  in_char$zip_ms$peak_finder$remove_double_peaks_in_scans()

  in_char$zip_ms$peak_finder$normalize_data()
  in_char$zip_ms$peak_finder$find_peaks_in_regions()
  in_char$zip_ms$peak_finder$add_offset()
  in_char$zip_ms$peak_finder$sort_ascending_mz()
  list(char_obj = in_char, processed = "doublenorm")
}

# remove points based on 99th percentile
# proper two pass normalization
# filter on frequency sd
filtersd = function(in_list){
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
  list(char_obj = in_char, processed = "filtersd")
}

write_peaks_for_assignment <- function(in_data){
  in_char = in_data$char_obj
  in_char$out_file = file.path(write_loc, paste0(in_data$processed, "_", in_char$zip_ms$peak_finder$sample_id, ".zip"))
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
    run_str = glue("pipenv run python3 ./Main.py /mlab/scratch/cesb_data/smirfe_dbs/n15_1600.db {zip_path} '_assigned.json' '[\"15N\"]' '[\"H\", \"Na\"]' '[1]'")
    system(run_str)
    setwd(curr_dir)
    out_value = return_file(assigned_file)
  } else {
    setwd("~/Projects/work/SMIRFE/SMIRFE_assigner/")
    run_str = glue("pipenv run python3 ./Main.py /mlab/scratch/cesb_data/smirfe_dbs/none_1600.db {zip_path} '_assigned.json' '[]' '[\"H\", \"Na\", \"K\", \"NH4\"]' '[1]'")
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
  sha256 = system2("sha256sum", args = in_file, stdout = TRUE)
  list(file = in_file, sha256 = sha256)
}

normalization_factors = function(...){
  processed_data = list(...)
  names(processed_data) = purrr::map_chr(processed_data, "processed")
  normed_samples = grep("nonorm", names(processed_data), value = TRUE, invert = TRUE)
  processed_data = processed_data[normed_samples]

  get_factors = function(processed_obj){
    norm_df = processed_obj$char_obj$zip_ms$peak_finder$peak_regions$normalization_factors
    norm_df$processing = processed_obj$processed
    norm_df
  }

  normalization_values = purrr::map_df(processed_data, get_factors)
  normalization_values

}

rsd_info = function(...) {
  processed_data = list(...)
  names(processed_data) = purrr::map_chr(processed_data, "processed")

  calc_rsd = function(in_char, processing){
    scan_peaks = purrr::map(in_char$zip_ms$peak_finder$peak_regions$peak_region_list, "peaks")
    n_peak = purrr::map_int(scan_peaks, nrow)
    keep_peak = n_peak >= 3
    scan_peaks = scan_peaks[keep_peak]
    rsd = purrr::map_df(scan_peaks, function(in_peaks){
      data.frame(mean = mean(in_peaks$Height),
                 sd = sd(in_peaks$Height),
                 n_peak = nrow(in_peaks),
                 stringsAsFactors = FALSE)
    })
    rsd$rsd = rsd$sd / rsd$mean
    rsd$processed = processing
    rsd
  }

  rsd_values = purrr::map_df(processed_data, function(in_data){
    calc_rsd(in_data$char_obj, in_data$processed)
  })

  rsd_values
}

summarize_rsd = function(rsd_df){
  get_mode = function(in_values){
    d_estimate = stats::density(in_values)
    y_max = which.max(d_estimate$y)
    d_estimate$x[y_max]
  }
  dplyr::group_by(rsd_df, processed) %>%
    dplyr::summarize(mean = mean(rsd),
                     median = median(rsd),
                     mode = get_mode(rsd),
                     max = max(rsd))
}

split_regions = function(p99_nonorm_data){
  # loadd("method_perc99_nonorm_data")
  # p99_nonorm_data = method_perc99_nonorm_data


  self = p99_nonorm_data$char_obj$zip_ms$peak_finder
  use_regions <- seq_len(length(self$peak_regions$peak_regions))
  signal_regions = self$peak_regions$peak_regions[use_regions]
  frequency_point_regions = self$peak_regions$frequency_point_regions
  tiled_regions = self$peak_regions$tiled_regions
  min_scan = self$peak_regions$min_scan
  peak_method = self$peak_method
  min_points = self$min_points

  signal_list = as.list(split(signal_regions, seq(1, length(signal_regions))))
  min_scan2 = floor(min_scan / 2)

  set.seed(1234)
  signal_list = signal_list[sample(length(signal_list), 500)]
  point_regions_list = purrr::map(signal_list, function(in_region){
    points = IRanges::subsetByOverlaps(frequency_point_regions, in_region)

    nonzero = as.data.frame(points@elementMetadata) %>%
      dplyr::filter(intensity > 0) %>%
      dplyr::pull(scan) %>% unique(.) %>% length(.)

    if (nonzero >= min_scan2) {
      tiles = IRanges::subsetByOverlaps(tiled_regions, in_region)
      return(list(points = points, tiles = tiles, region = in_region))
    } else {
      return(NULL)
    }
  })

  not_null = purrr::map_lgl(point_regions_list, ~ !is.null(.x))
  point_regions_list = point_regions_list[not_null]
  point_regions_list = point_regions_list[sample(length(point_regions_list))]

  # split_data = purrr::imap(point_regions_list, function(.x, .y){
  #   message(.y)
  #   split_region_by_peaks(.x, peak_method = peak_method, min_points = min_points)
  # })
  split_data = FTMS.peakCharacterization:::internal_map$map_function(point_regions_list, split_region_by_peaks,
                                         peak_method = peak_method, min_points = min_points)

  null_regions = purrr::map_lgl(split_data, ~ is.null(.x[[1]]$points))
  split_data = split_data[!null_regions]


  n_regions = purrr::map_int(split_data, ~length(.x))

  possible_region = which(n_regions == 2)
  possible_region = which(n_regions == 2)
  n_each = purrr::imap_dfr(split_data[possible_region], function(.x, .y){
    data.frame(region = .y, n_1 = nrow(.x[[1]]$peaks), n_2 = nrow(.x[[2]]$peaks))
  })
  examine_region = dplyr::filter(n_each, n_1 > 50, n_2 > 50)
  use_list = point_regions_list[[examine_region[1, "region"]]]

  use_list

}

density_calculate = function(frequency_values, metadata, frequency_range = NULL){
  freq_points = FTMS.peakCharacterization::frequency_points_to_frequency_regions(as.data.frame(frequency_values), multiplier = metadata$multiplier)

  sliding_regions = FTMS.peakCharacterization::create_frequency_regions(point_spacing = metadata$point_spacing, frequency_range = frequency_range,
                                             n_point = metadata$n_point, delta_point = metadata$delta_point,
                                             multiplier = metadata$multiplier)

  sliding_density = IRanges::countOverlaps(sliding_regions, freq_points)

  n_region = 9 * metadata$n_point
  small_window = c(1 * metadata$n_point, 5 * metadata$n_point)
  large_window = c((5 * metadata$n_point) + 1, 9 * metadata$n_point)

  start_i = (9 * metadata$n_point) + 1
  end_i = length(sliding_regions) - (9 * metadata$n_point) - 1

  stats = purrr::map_dbl(seq(start_i, end_i), function(in_i){
    inner_window = list(sort(in_i - small_window),
                        sort(in_i + small_window))
    outer_window = list(sort(in_i - large_window),
                        sort(in_i + large_window))

    inner_seq = unlist(purrr::map(inner_window, ~ seq(.x[1], .x[2])))
    outer_seq = unlist(purrr::map(outer_window, ~ seq(.x[1], .x[2])))

    calculate_statistic = function(query_density, reference_window){
      ref_mean = mean(reference_window)
      ref_sd = sd(reference_window)
      if ((ref_mean > 0) & (ref_sd > 0)) {
        z = ((query_density - ref_mean)^2) / (ref_sd^2)
        return(z * sign(query_density - ref_mean))
      } else {
        return(0)
      }

    }

    z_inner = calculate_statistic(sliding_density[in_i], sliding_density[inner_seq])
    z_outer = calculate_statistic(sliding_density[in_i], sliding_density[outer_seq])
    max(z_inner, z_outer)
  })
  all_z = rep(0, length(sliding_density))
  all_z[seq(start_i, end_i)] = stats

  tmp_df = data.frame(z = all_z, frequency_start = IRanges::start(sliding_regions))

  high_df = which(tmp_df$z >= 100)

  high_sliding = sliding_regions[high_df]
  high_merged = IRanges::reduce(high_sliding)
  min_width = IRanges::width(sliding_regions[1]) / 10 * 3
  keep_merged = high_merged[IRanges::width(high_merged) >= min_width]
  hpd_points = purrr::map_df(seq(1, length(keep_merged)), function(in_region){
    tmp_points = IRanges::subsetByOverlaps(freq_points, keep_merged[in_region])
    tmp_data = as.data.frame(tmp_points@elementMetadata)
    tmp_data$region = in_region
    tmp_data
  })
  return(list(hpd_points = hpd_points, freq_points = freq_points, hpd_regions = keep_merged, stats = tmp_df))
}

hpds_from_excel = function(in_data){
  sample = in_data$char_obj$zip_ms$id
  match_excel = grep(sample, excel_files, value = TRUE)

  xl_data = readxl::read_excel(match_excel, skip = 8, col_names = FALSE)
  xl_data = xl_data[, 1:2]
  names(xl_data) = c("mz", "intensity")
  message("got xl data")
  frequency_coefficients = in_data$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions@metadata$frequency_coefficients
  frequency_description = in_data$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions@metadata$frequency_fit_description
  mz_coefficients = in_data$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions@metadata$mz_coefficients
  mz_description = in_data$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions@metadata$mz_fit_description
  org_frequency = as.data.frame(in_data$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions@elementMetadata)
  xl_data$frequency = FTMS.peakCharacterization:::predict_exponentials(xl_data$mz, frequency_coefficients, frequency_description)

  sliding_metadata = in_data$char_obj$zip_ms$peak_finder$peak_regions$sliding_regions@metadata

  freq_diff = 1000
  sliding_metadata$point_spacing = round(freq_diff / 10)


  frequency_range = c(min(c(min(c(xl_data$frequency, org_frequency$frequency)))),
                      max(c(max(c(xl_data$frequency, org_frequency$frequency)))))

  xl_hpd = density_calculate(xl_data, sliding_metadata, frequency_range)

  peak_data = in_data$char_obj$zip_ms$peak_finder$peak_regions$peak_data
  scan_level = in_data$char_obj$zip_ms$peak_finder$peak_regions$scan_level_arrays
  split_xl = split(xl_hpd$hpd_points, xl_hpd$hpd_points$region)
  xl_ranges = purrr::map(split_xl, ~ range(.x$mz))
  message("got hpd sites")
  message(nrow(peak_data))
  inside_peaks = purrr::imap_dfr(xl_ranges, function(in_range, in_id){
    tmp_df = dplyr::filter(peak_data, dplyr::between(ObservedMZ, in_range[1], in_range[2]))
    if (nrow(tmp_df) > 0) {
      tmp_df$region = in_id
      return(tmp_df)
    } else {
      return(NULL)
    }
  })
  message("inside")
  list(sample = sample, peak_data = inside_peaks, hpd = xl_hpd,
       processed = in_data$processed)
}

