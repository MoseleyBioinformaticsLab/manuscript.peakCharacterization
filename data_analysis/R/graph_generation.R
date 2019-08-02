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

plot_frequency_conversion = function(raw_data){
  raw_points = as.data.frame(raw_data$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions@elementMetadata)
  raw_scan = dplyr::filter(raw_points, scan %in% unique(scan)[1])
  raw_points2 = convert_mz_frequency(raw_scan)
  convertable_stretch = rle(raw_points2$convertable)
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
  raw_possible = purrr::map(possible_index, ~ raw_points2[.x, ])
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

  all_frequency_difference_plot = ggplot(raw_points2, aes(x = mz, y = log10(abs(mean_freq_diff)))) + geom_point() +
    geom_point(data = dplyr::filter(raw_points2, convertable), color = "red") + labs(x = "M/Z", y = "Log10(Frequency Differences)")

  convertable_raw2 = dplyr::filter(raw_points2, convertable)
  raw2_model = FTMS.peakCharacterization:::fit_exponentials(convertable_raw2$mean_mz, convertable_raw2$mean_frequency, c(0, -1/2, -1/3))
  convertable_raw2$predict_frequency = FTMS.peakCharacterization:::predict_exponentials(convertable_raw2$mean_mz, raw2_model$coefficients, c(0, -1/2, -1/3))
  mz_frequency_plot = ggplot(convertable_raw2, aes(x = mean_mz, y = mean_frequency)) + geom_point(size = 4) +
    geom_line(aes(y = predict_frequency), color = "red", size = 2) +
    labs(x = "M/Z", y = "Frequency")

  sqrt_model = FTMS.peakCharacterization:::fit_exponentials(convertable_raw2$mean_mz, convertable_raw2$mean_frequency, c(0, -1/2))
  sqrt_pred = FTMS.peakCharacterization:::predict_exponentials(convertable_raw2$mean_mz, sqrt_model$coefficients, c(0, -1/2))

  residual_data = rbind(data.frame(mz = convertable_raw2$mean_mz, frequency = convertable_raw2$mean_frequency,
                                   predicted = convertable_raw2$predict_frequency,
                                   type = "cube", stringsAsFactors = FALSE),
                        data.frame(mz = convertable_raw2$mean_mz, frequency = convertable_raw2$mean_frequency,
                                   predicted = sqrt_pred,
                                   type = "square", stringsAsFactors = FALSE)) %>%
    dplyr::mutate(residuals = (predicted - frequency) / frequency)
  compare_models = ggplot(residual_data, aes(x = mz, y = residuals)) + geom_point() +
    facet_wrap(~ type, ncol = 1)

  frequency_plot =  ((mz_point_plot / frequency_point_plot / point_difference_plot)  |
    (all_frequency_difference_plot / mz_frequency_plot)) + plot_annotation(tag_levels = "A")
  frequency_plot

}
