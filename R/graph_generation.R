plot_nap_peaks = function(raw_data, processed_data, assignments, interesting_regions){
  sample_ids = c(raw_data$char_obj$peak_finder$sample_id,
                 processed_data$peak_finder$sample_id,
                 assignments$sample)
  if (length(unique(sample_ids)) > 1) {
    stop("There are multiple samples in this data collection!")
  }

  region1 = interesting_regions$emfs[[1]]$peak_info[[1]]
  region1$adduct_emf = paste0(region1$adduct_IMF, region1$complete_EMF)
  region_split = split(region1, region1$adduct_emf)[[1]]

  region_peaks = unique(region_split$PeakID)
  assignment_data = dplyr::filter(assignments$data, PeakID %in% region_peaks)

  region_split = tidyr::spread(region_split, key = Type, value = Assignment_Data)

  region_split = dplyr::mutate(region_split, NAP = as.numeric(NAP), relnap = NAP / max(NAP))

  assignment_data = tidyr::spread(assignment_data, key = Measurement, value = Value)

  region_split = dplyr::left_join(region_split, assignment_data, by = "PeakID")
  region_split = dplyr::mutate(region_split, relIntensity = log10(relnap * max(Height)))

  peak_regions = processed_data$zip_ms$peak_finder$peak_regions

  data_regions = dplyr::filter(peak_regions$peak_data, PeakID %in% region_peaks)
  scan_regions = peak_regions$scan_peaks[data_regions$PeakID]
  scan_data = purrr::imap_dfr(scan_regions, function(.x, .y){
    .x$peak = .y
    as.data.frame(.x)
  })

  scan_means = dplyr::group_by(scan_data, peak) %>%
    dplyr::summarise(mean_mz = mean(ObservedMZ), mean_height = mean(RawHeight))

  xcal_file = file.path("data_analysis/data_input", paste0(sample_ids[1], ".xlsx"))
  xcal_data = readxl::read_xlsx(xcal_file, skip = 6) %>% dplyr::mutate(mz = m/z)

  xcal_match = purrr::map_df(seq(1, nrow(region_split)), function(in_row){
    min_loc = which.min(abs(xcal_data$mz - region_split[in_row, "ObservedMZ"]))
    xcal_data[min_loc, ]
  })

  ggplot(scan_data, aes(x = ObservedMZ, y = log10(RawHeight))) + geom_point() +
    geom_point(data = region_split, aes(x = ObservedMZ, y = relIntensity), color = "red", size = 4) +
    geom_point(data = xcal_match, aes(x = mz-0.02, y = log10(Intensity)), color = "lightblue", size = 4) +
    geom_point(data = scan_means, aes(x = mean_mz + 0.02, y = log10(mean_height)), color = "brown", size = 4)
}

plot_frequency_conversion = function(char_obj){
  raw_points = as.data.frame(char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions$frequency[[1]]@elementMetadata)
  is_convertable = dplyr::between(raw_points$mean_freq_diff, 0.49, 0.51)
  is_convertable[is.na(is_convertable)] = FALSE
  raw_points$org_convertable = is_convertable
  convertable_stretch = rle(raw_points$org_convertable)
  convertable_df = data.frame(values = convertable_stretch$values,
                              lengths = convertable_stretch$lengths) %>%
    dplyr::mutate(index = seq(1, length(convertable_stretch$lengths)))

  index_list = vector("list", length = nrow(convertable_df))
  start_index = 1
  for (i in seq_along(index_list)) {
    end_index = start_index + convertable_df$lengths[i] - 1
    index_list[[i]] = seq(start_index, end_index)
    start_index = end_index + 1
  }

  possible_index = index_list[dplyr::filter(convertable_df, values, lengths > 4) %>% dplyr::pull(index)]
  raw_possible = purrr::map(possible_index, ~ raw_points[.x, ])
  raw_high = purrr::map_lgl(raw_possible, ~ sum(.x$intensity >= 1e4) > 1)

  example_data = raw_possible[raw_high][[1]]

  example_data$pair_mz = example_data$mean_mz
  example_data$pair_difference = example_data$mean_offset

  example_data$pair_frequency = example_data$mean_frequency
  example_data$pair_frequency_offset = example_data$mean_freq_diff

  example_data$level1_intensity = -.1
  example_data$level1_group = 0.1
  example_data$mz_group_xstart = example_data$mz + 1e-5
  example_data$mz_group_xend = dplyr::lead(example_data$mz - 1e-5)
  example_data$frequency_group_xstart = example_data$mean_frequency - 5e-2
  example_data$frequency_group_xend = dplyr::lead(example_data$mean_frequency) + 5e-2
  example_data[1, "pair_mz"] = NA
  mz_point_plot = ggplot(example_data, aes(x = mz, y = log10(intensity+1))) + geom_point() +
    geom_segment(aes(x = mz_group_xstart, xend = mz_group_xend, y = level1_intensity, yend = level1_intensity), color = "red") +
    geom_point(aes(x = pair_mz, y = level1_intensity), color = "red") + labs(y = "Log10(Intensity)", x = "M/Z")

  frequency_point_plot = ggplot(example_data, aes(x = pair_frequency, y = log10(intensity+1))) + geom_point() +
    labs(y = "Log10(Intensity)", x = "Frequency") + geom_segment(aes(x = frequency_group_xstart, xend = frequency_group_xend, y = level1_intensity, yend = level1_intensity), color = "red")

  point_difference_plot = ggplot(example_data, aes(x = pair_frequency, y = pair_frequency_offset)) + geom_point() +
    labs(x = "Frequency", y = "Frequency Differences")

  all_frequency_difference_plot = ggplot(raw_points, aes(x = mz, y = log10(abs(mean_freq_diff)))) + geom_point() +
    geom_point(data = dplyr::filter(raw_points, org_convertable), color = "red") + labs(x = "M/Z", y = "Log10(Frequency Differences)")

  convertable_raw2 = dplyr::filter(raw_points, org_convertable)
  raw2_model = FTMS.peakCharacterization:::fit_exponentials(convertable_raw2$mean_mz, convertable_raw2$mean_frequency, c(0, -1/2, -1/3))
  convertable_raw2$predict_frequency = FTMS.peakCharacterization:::predict_exponentials(convertable_raw2$mean_mz, raw2_model$coefficients, c(0, -1/2, -1/3))
  mz_frequency_plot = ggplot(convertable_raw2, aes(x = mean_mz, y = mean_frequency)) + geom_point(size = 4) +
    geom_line(aes(y = predict_frequency), color = "red", size = 2) +
    labs(x = "M/Z", y = "Frequency")



  sqrt_model = FTMS.peakCharacterization:::fit_exponentials(convertable_raw2$mean_mz, convertable_raw2$mean_frequency, c(0, -1/2))
  sqrt_pred = FTMS.peakCharacterization:::predict_exponentials(convertable_raw2$mean_mz, sqrt_model$coefficients, c(0, -1/2))

  convertable_raw2$predict_square = sqrt_pred
  sqrt_frequency_plot = ggplot(convertable_raw2, aes(x = mean_mz, y = mean_frequency)) + geom_point(size = 4) +
    geom_line(aes(y = predict_square), color = "red", size = 2) +
    labs(x = "M/Z", y = "Frequency")
  residual_data = rbind(data.frame(mz = convertable_raw2$mean_mz, frequency = convertable_raw2$mean_frequency,
                                   predicted = convertable_raw2$predict_frequency,
                                   type = "cube", stringsAsFactors = FALSE),
                        data.frame(mz = convertable_raw2$mean_mz, frequency = convertable_raw2$mean_frequency,
                                   predicted = sqrt_pred,
                                   type = "square", stringsAsFactors = FALSE)) %>%
    dplyr::mutate(residuals = (predicted - frequency) / frequency)

  frequency_plot =  ((mz_point_plot / frequency_point_plot / point_difference_plot)  |
    (all_frequency_difference_plot / mz_frequency_plot)) + plot_annotation(tag_levels = "A")
  frequency_histogram = raw_points %>%
    ggplot(aes(x = log10(mean_freq_diff))) +
    geom_histogram(bins = 100) +
    labs(x = "Log10(Frequency Differences)")
  list(all = frequency_plot, histrogram = frequency_histogram,
       residuals = residual_data)

}

plot_peak_ordering = function(data_obj){
  peak_data = data_obj$char_obj
  scan_peaks = purrr::map(peak_data$zip_ms$peak_finder$peak_regions$peak_region_list, "peaks")
  scan_level_peaks = as.data.frame(scan_peaks[[1]])

  scan_models = peak_data$zip_ms$peak_finder$peak_regions$frequency_point_regions$metadata$frequency_coefficients_all
  scan_models = dplyr::filter(scan_models, scan %in% scan_level_peaks$scan)

  by_scan = purrr::map_df(seq(1, nrow(scan_level_peaks)), function(in_row){
    use_mz = scan_level_peaks[in_row, "ObservedMZ"]
    use_model = dplyr::filter(scan_models, scan %in% scan_level_peaks[in_row, "scan"])

    out_freq = FTMS.peakCharacterization:::predict_exponentials(use_mz, unlist(use_model[1, c("V1", "V2", "V3")]),
                                                                peak_data$zip_ms$peak_finder$peak_regions$frequency_point_regions$metadata$frequency_fit_description)
    data.frame(mz = use_mz, single_frequency = scan_level_peaks[in_row, "ObservedFrequency"],
               scan_frequency = out_freq)
  })

  by_scan = dplyr::arrange(by_scan, mz) %>% dplyr::mutate(mz_order = seq(1, nrow(.))) %>%
    dplyr::arrange(single_frequency) %>% dplyr::mutate(single_order = seq(1, nrow(.))) %>%
    dplyr::arrange(scan_frequency) %>% dplyr::mutate(scan_order = seq(1, nrow(.)))

  by_scan = tidyr::gather(by_scan, key = "type_order", value = "order", single_order, scan_order)

  peak_ordering_plot = ggplot(by_scan, aes(x = mz_order, y = order, color = type_order)) + geom_point(size = 3) +
    theme(legend.position = c(0.8, 0.9))
  peak_ordering_plot
}

plot_sliding_window_density = function(peak_data){
  sliding_regions = peak_data$zip_ms$peak_finder$peak_regions$sliding_regions
  point_region_list = peak_data$zip_ms$peak_finder$peak_regions$frequency_point_regions$frequency
  multiplier = peak_data$zip_ms$peak_finder$quantile_multiplier
  n_point_region = peak_data$zip_ms$peak_finder$n_point_region

  nz_counts = FTMS.peakCharacterization:::count_overlaps(sliding_regions, point_region_list[[1]])
  n_region = seq(2, length(point_region_list))
  for (iregion in n_region) {
    nz_counts_iregion = FTMS.peakCharacterization:::count_overlaps(sliding_regions, point_region_list[[iregion]])
    nz_counts = nz_counts + nz_counts_iregion
  }

  chunk_indices = seq(1, length(nz_counts), by = n_point_region)

  chunk_perc = purrr::map_dbl(chunk_indices, function(in_index){
    use_counts = nz_counts[seq(in_index, min(in_index + (n_point_region - 1), length(nz_counts)))]
    if (max(use_counts) > 0) {
      return(stats::quantile(use_counts, 0.99))
    } else {
      return(0)
    }
  })

  single_chunk = data.frame(use_counts = nz_counts[
    seq(chunk_indices[1], (chunk_indices[2] - 1))])
  single_99 = stats::quantile(single_chunk$use_counts, 0.99)

  all_99 = data.frame(perc_99 = chunk_perc,
                           index = seq(1, length(chunk_perc)))
  all_cut = 1.5 * median(all_99$perc_99)
  single_plot = ggplot(single_chunk, aes(x = use_counts)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = single_99, color = "blue") +
    scale_y_log10(expand = c(0, 0), limits = c(1, NA)) +
    labs(x = "Number of non-zero points",
         y = "Log10(counts)")
  all_plot = ggplot(all_99, aes(x = perc_99)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = all_cut, color = "red") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "99th percentiles across segments")
  (single_plot | all_plot) + plot_annotation(tag_levels = "A")
}

plot_rsd_differences = function(rsd_values){

  rsd_values$processed = forcats::fct_relevel(rsd_values$processed, "noperc_nonorm", "perc99_nonorm", "singlenorm", "singlenorm_int", "doublenorm",  "filtersd")
  rsd_comparison_plot = ggplot(rsd_values, aes(x = rsd, y = processed, fill = processed)) + geom_density_ridges() +
    coord_cartesian(xlim = c(0, 1)) + theme(legend.position = "none")
  rsd_comparison_plot
}

proportional_error = function(peak_data){
  raw_peaks = peak_data$peak_finder$peak_regions$scan_peaks
  n_peak = purrr::map_int(raw_peaks, nrow)
  raw_peaks = raw_peaks[n_peak > 2]

  peak_sd = purrr::map_df(raw_peaks, function(in_peaks){
    tmp_df = data.frame(mean = mean(in_peaks$Height),
                        sd = sd(in_peaks$Height),
                        stringsAsFactors = FALSE)
  })
  peak_sd$rsd = peak_sd$sd / peak_sd$mean
}


normalization_graph = function(normalization_data){
  normalization_data = dplyr::filter(normalization_data, !(processing %in% "filtersd"))

  hist_plot = ggplot(normalization_data, aes(x = normalization)) + geom_histogram() +
    facet_wrap(~ processing, ncol = 1)

  split_norm = split(normalization_data, normalization_data$processing)
  ref_norm = split_norm[["singlenorm"]]
  diff_norm = purrr::map_df(split_norm, function(in_norm){
    combine_norm = dplyr::left_join(ref_norm, in_norm, by = "scan", suffix = c(".ref", ".in"))
    dplyr::mutate(combine_norm, scan = scan, processing = processing.in,
                  diff = normalization.in - normalization.ref) %>%
      dplyr::select(scan, diff, processing)
  })

  diff_norm = dplyr::filter(diff_norm, !(processing %in% "singlenorm"))
  diff_plot = ggplot(diff_norm, aes(x = scan, y = diff, color = processing)) + geom_point() + geom_line() +
    theme(legend.position = c(0.4, 0.9))

  (hist_plot | diff_plot) + plot_annotation(tag_levels = "A")
}

correlate_scan_height = function(...){
  processed_data = list(...)
  names(processed_data) = purrr::map_chr(processed_data, "processed")
  use_names = grep("perc99", names(processed_data), value = TRUE)
  message(use_names)
  processed_obj = processed_data[[use_names]]
  # to make this something we can test
  # loadd("method_perc99_nonorm_data")
  # processed_obj = method_perc99_nonorm_data

  intensity_measure = c("RawHeight", "Height")
  summary_function = median
  normalize_peaks = "both"

  scan_peaks = purrr::map(processed_obj$char_obj$zip_ms$peak_finder$peak_regions$peak_region_list, ~ as.data.frame(.x$peaks))

  ## getting the relationship between normalized peak height and scan ----
  normalization_factors <- FTMS.peakCharacterization:::single_pass_normalization(scan_peaks, intensity_measure = intensity_measure, summary_function = summary_function,
                                                                                 min_ratio = 0)

  normed_peaks = purrr::map(scan_peaks, function(.x){
    .x2 = dplyr::inner_join(.x, normalization_factors, by = "scan")
    .x2$Height = exp(log(.x2$RawHeight) - .x2$normalization)
    .x2$normalization = NULL
    .x2
  })

  peak_scan_correlation = purrr::imap_dfr(normed_peaks, function(.x, .y){
    data.frame(n_scan = nrow(.x), correlation = cor(log10(.x$Height), .x$scan),
               peak = as.numeric(.y), stringsAsFactors = FALSE)
  }) %>% dplyr::filter(!is.na(correlation))

  max_scan = max(dplyr::filter(peak_scan_correlation, n_scan <= 150)$n_scan)
  high_cor = dplyr::filter(peak_scan_correlation, n_scan >= 0.8*max_scan, correlation >= 0.8)

  example_peak = normed_peaks[[high_cor[1, "peak"]]]

  height_2_scan = ggplot(example_peak, aes(x = scan, y = log10(Height))) + geom_point(size = 3)


  # do the basic normalization to get M/Z to difference ----
  intensity_measure = c("RawHeight", "Height")
  summary_function = median
  use_peaks = NULL
  min_ratio = 0
  n_scan_per_peak <- purrr::map_int(scan_peaks, function(x){
    if (sum(duplicated(x$scan)) == 0) {
      return(length(x$scan))
    } else {
      return(0L)
    }
  })

  use_measure <- NULL
  for (imeasure in intensity_measure){
    if (imeasure %in% names(scan_peaks[[1]])) {
      use_measure <- imeasure
      break()
    }
  }

  if (is.null(use_measure)) {
    stop("The desired intensity measure for normalization is not present!")
  }

  scan_cutoff <- quantile(n_scan_per_peak, 0.95)

  if (is.null(use_peaks)) {
    use_peaks <- rep(TRUE, length(scan_peaks))
  }

  normalize_peaks <- which((n_scan_per_peak >= scan_cutoff) & use_peaks)

  all_scans <- data.frame(scan = unique(unlist(purrr::map(scan_peaks, function(x){unique(x$scan)}))))

  peak_intensity <- purrr::map_dfc(scan_peaks[normalize_peaks], function(x){
    tmp_data <- dplyr::left_join(all_scans, as.data.frame(x[, c("scan", use_measure)]), by = "scan")
    log(tmp_data[, use_measure])
  })

  peak_mz = purrr::map_dfc(scan_peaks[normalize_peaks], function(x){
    tmp_data <- dplyr::left_join(all_scans, as.data.frame(x[, c("scan", "ObservedMZ")]), by = "scan")
    tmp_data[, "ObservedMZ"]
  })
  all_na = purrr::map_lgl(seq_len(nrow(peak_intensity)), ~ sum(is.na(peak_intensity[.x, ])) == ncol(peak_intensity))

  all_scans = all_scans[!all_na, , drop = FALSE]
  peak_intensity = peak_intensity[!all_na, , drop = FALSE]
  peak_mz = peak_mz[!all_na, , drop = FALSE]

  intensity_ratio <- purrr::map_dfr(seq_len(nrow(peak_intensity)), function(x){
    #message(x)
    peak_intensity[x, ] / max(peak_intensity[x, ], na.rm = TRUE)
  })
  peak_intensity[intensity_ratio < min_ratio] <- NA

  intensity_scans <- purrr::map_int(seq_len(ncol(peak_intensity)), function(x){
    sum(!is.na(peak_intensity[, x]))
  })

  peak_intensity <- peak_intensity[, intensity_scans >= scan_cutoff, drop = FALSE]
  peak_mz = peak_mz[, intensity_scans >= scan_cutoff, drop = FALSE]

  notna_scans <- rowSums(!is.na(as.matrix(peak_intensity)))
  keep_scans <- notna_scans >= 25

  if (sum(keep_scans) == 0) {
    stop("No scans left to use in normalization!")
  }

  all_scans <- all_scans[keep_scans, , drop = FALSE]
  peak_intensity <- peak_intensity[keep_scans, ]
  peak_mz = peak_mz[keep_scans, ]

  scan_distances <- purrr::map_dbl(seq(1, nrow(peak_intensity)), function(in_scan){
    scan_peaks <- peak_intensity[in_scan, ,drop = FALSE]
    scan_peaks_matrix <- matrix(unlist(scan_peaks), nrow = nrow(peak_intensity) - 1, ncol = ncol(scan_peaks), byrow = TRUE)

    other_matrix <- as.matrix(peak_intensity[-in_scan, , drop = FALSE])
    scan_other_diff <- scan_peaks_matrix - other_matrix
    scan_distances <- purrr::map_dbl(seq(1, nrow(scan_other_diff)), function(x){
      sqrt(sum(scan_other_diff[x, ]^2, na.rm = TRUE))
    })
    sum(scan_distances)
  })

  normalize_scan <- which.min(scan_distances)

  scan_norm_matrix <- matrix(unlist(peak_intensity[normalize_scan, , drop = FALSE]),
                             nrow = nrow(peak_intensity), ncol = ncol(peak_intensity), byrow = TRUE)

  diff_matrix <- as.matrix(peak_intensity) - scan_norm_matrix

  rownames(diff_matrix) = paste0("d", seq(1, nrow(diff_matrix)))
  colnames(diff_matrix) = paste0("p", seq(1, ncol(diff_matrix)))

  diff_df = as.data.frame(t(diff_matrix))
  diff_df$peak = rownames(diff_df)
  diff_df$mz = colMeans(peak_mz, na.rm = TRUE)

  int_df = as.data.frame(t(as.matrix(peak_intensity)))
  colnames(int_df) = paste0("i", seq(1, ncol(int_df)))
  int_df$peak = diff_df$peak

  diff_int_df = dplyr::left_join(diff_df, int_df, by = "peak")

  diff_int_list = purrr::map(seq(1, nrow(diff_matrix)), function(in_scan){
    grab_cols = paste0(c("d", "i"), in_scan)
    tmp_df = diff_int_df[, grab_cols]
    names(tmp_df) = c("difference", "intensity")
    tmp_df = dplyr::filter(tmp_df, !is.na(difference), !is.na(intensity))
    tmp_df$use = FALSE
    cutoff = max(tmp_df$intensity) * 0.7
    tmp_df[tmp_df$intensity >= cutoff, "use"] = TRUE
    tmp_df
  })

  diff_int_plots = purrr::map(diff_int_list, function(intensity_df){
    p = ggplot(dplyr::filter(intensity_df, !use), aes(x = intensity, y = difference)) + geom_point(size = 3) +
      geom_point(data = dplyr::filter(intensity_df, use), size = 3, color = "red")
    p
  })

  list(correlation = height_2_scan, diff_df = diff_int_list, diff_plot = diff_int_plots)

}

correlate_scan_height_graph = function(correlation_plots){
  p = (correlation_plots$correlation | correlation_plots$diff_plot[[1]]) + plot_annotation(tag_levels = "A")
  p
}


plot_peak_fitting = function(processed_obj){
  # to make this something we can test
  #loadd("method_perc99_nonorm_data")
  #processed_obj = method_perc99_nonorm_data

  use_peak = processed_obj$char_obj$zip_ms$peak_finder$peak_regions$peak_region_list[[1]]
  peak_points = as.data.frame(use_peak$points@elementMetadata)
  use_scan = unique(peak_points$scan)[1]
  peak_points = dplyr::filter(peak_points, scan %in% use_scan)

  w = peak_points$intensity / max(peak_points$intensity)
  peak_points$log_int = log(peak_points$intensity + 1e-8)
  peak_fit = FTMS.peakCharacterization:::parabolic_fit(peak_points$frequency, peak_points$log_int, w = w)
  peak_points$log_fit = peak_fit$fitted.values
  fitted_info = FTMS.peakCharacterization:::get_fitted_peak_info(peak_points, "frequency", w)

  log_plot = ggplot(peak_points, aes(x = frequency, y = log_int)) +
    geom_line(aes(x = frequency, y = log_fit), color = "red", size = 2) +
    geom_point(size = 3) +
    geom_point(data = fitted_info, aes(x = ObservedCenter, y = log(Height)), color = "blue", size = 3) +
    labs(x = "Frequency", y = "Log(Intensity)")

  norm_plot = ggplot(peak_points, aes(x = frequency, y = intensity)) +
    geom_line(aes(x = frequency, y = exp(log_fit)), color = "red", size = 2) +
    geom_point(size = 3) +
    geom_point(data = fitted_info, aes(x = ObservedCenter, y = Height), color = "blue", size = 3) +
    labs(x = "Frequency", y = "Intensity")

  (log_plot | norm_plot) + plot_annotation(tag_levels = "A")
}

plot_region_splitting = function(region_list){
  frequency_point_regions = region_list$points
  tiled_regions = region_list$tiles
  frequency_point_regions@elementMetadata$log_int <- log(frequency_point_regions@elementMetadata$intensity + 1e-8)

  scan_runs = rle(frequency_point_regions@elementMetadata$scan)
  #if (min(scan_runs$lengths) < min_points + 2) {
  frequency_point_splitscan <- split(frequency_point_regions, frequency_point_regions@elementMetadata$scan)
  reduced_peaks <- purrr::map_df(names(frequency_point_splitscan), function(in_scan){
    FTMS.peakCharacterization:::get_reduced_peaks(frequency_point_splitscan[[in_scan]], peak_method = "lm_weighted", min_points = 4)
  })

  frequency_multiplier = region_list$points@metadata$frequency_multiplier
  reduced_regions = IRanges::IRanges(start = round(reduced_peaks$ObservedCenter.frequency * frequency_multiplier),
                                     width = 1)
  tiled_overlap = IRanges::countOverlaps(tiled_regions, reduced_regions)

  point_data = as.data.frame(region_list$points@elementMetadata)


  tiled_regions = as.data.frame(region_list$tiles)
  tiled_regions$start_freq = tiled_regions$start / frequency_multiplier
  tiled_regions$end_freq = tiled_regions$end / frequency_multiplier
  tiled_regions$counts = tiled_overlap
  tiled_regions$mid_point = tiled_regions$start_freq + 0.25
  tiled_regions$intensity = 10000

  original_points = as.data.frame(frequency_point_regions@elementMetadata)

  xmin = min(original_points$frequency - 0.5)
  xmax = max(original_points$frequency + 0.5)

  p1 = ggplot(original_points, aes(x = frequency, y = intensity, group = scan)) +
    geom_point() + geom_line() + coord_cartesian(xlim = c(xmin, xmax)) +
    labs(x = NULL)

  p2 = ggplot(reduced_peaks, aes(x = ObservedCenter.frequency, y = Height.frequency)) +
    geom_point() + geom_segment(data = tiled_regions, aes(x = start_freq + 5e-2, xend = end_freq - 5e-2, y = intensity, yend = intensity), color = "red") +
    coord_cartesian(xlim = c(xmin, xmax)) +
    labs(x = NULL, y = "intensity")

  p3 = ggplot(tiled_regions, aes(x = mid_point, y = counts)) + geom_col() +
    coord_cartesian(xlim = c(xmin, xmax)) +
    labs(x = "frequency", y = "no. of peaks")

  (p1 / p2 / p3) + plot_annotation(tag_level = "A")
}
