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


hpds_from_excel = function(in_data, excel_files){
  sample = in_data$char_obj$zip_ms$id
  match_excel = grep(sample, excel_files, value = TRUE)

  xl_data = readxl::read_excel(match_excel, skip = 8, col_names = FALSE)
  xl_data = xl_data[, 1:2]
  names(xl_data) = c("mz", "intensity")
  message("got xl data")
  frequency_coefficients = in_data$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions$metadata$frequency_coefficients
  frequency_description = in_data$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions$metadata$frequency_fit_description
  mz_coefficients = in_data$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions$metadata$mz_coefficients
  mz_description = in_data$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions$metadata$mz_fit_description
  org_frequency_peaks = purrr::map_df(seq_len(length(in_data$char_obj$zip_ms$peak_finder$peak_regions$peak_region_list)),
                                      function(.x){
                                        tmp = as.data.frame(in_data$char_obj$zip_ms$peak_finder$peak_regions$peak_region_list[[.x]]$peaks[, c("ObservedMZ", "ObservedFrequency", "scan")])
                                        tmp$region = .x
                                        tmp
                                      })
  xl_data$frequency = FTMS.peakCharacterization:::predict_exponentials(xl_data$mz, frequency_coefficients, frequency_description)

  sliding_metadata = in_data$char_obj$zip_ms$peak_finder$peak_regions$sliding_regions@metadata

  freq_diff = 1000
  sliding_metadata$point_spacing = round(freq_diff / 10)


  frequency_range = c(min(c(min(c(xl_data$frequency, org_frequency_peaks$ObservedFrequency)))),
                      max(c(max(c(xl_data$frequency, org_frequency_peaks$ObservedFrequency)))))

  xl_hpd = density_calculate(xl_data, sliding_metadata, frequency_range)

  peak_data = in_data$char_obj$zip_ms$peak_finder$peak_regions$peak_data
  scan_level = in_data$char_obj$zip_ms$peak_finder$peak_regions$scan_level_arrays
  split_xl = split(xl_hpd$hpd_points, xl_hpd$hpd_points$region)
  xl_ranges = purrr::map(split_xl, ~ range(.x$frequency))
  message("got hpd sites")
  message(nrow(peak_data))
  inside_peaks = purrr::imap_dfr(xl_ranges, function(in_range, range_id){
    tmp_df = dplyr::filter(peak_data, dplyr::between(ObservedFrequency, in_range[1], in_range[2]))
    if (nrow(tmp_df) > 0) {
      tmp_df$HPD = TRUE
      tmp_df$HPDID = range_id
      return(tmp_df[, c("PeakID", "HPD", "HPDID")])
    } else {
      return(NULL)
    }
  })
  peak_data2 = dplyr::left_join(peak_data, inside_peaks, by = "PeakID")
  peak_data2[is.na(peak_data2$HPD), "HPD"] = FALSE
  #message("inside")
  max_sd = max(peak_data2$ObservedFrequencySD)
  list(peak_data = peak_data2, xl_data = xl_data, hpd = xl_hpd, hpd_ranges = xl_ranges,
       sample_id = rename_samples(sample))
}


chisq_hpds = function(hpd_res){
  peak_data = hpd_res$peak_data
  cont_table = table(peak_data[, c("HighSD", "HPD")])
  chisq_res = broom::tidy(chisq.test(cont_table))
  chisq_res$sample = hpd_res$sample
  chisq_res
}

width_sd_hpds = function(hpd_res){
  peak_data = hpd_res$peak_data
  split_hpd = split(peak_data, peak_data$HPDID)
  tmp_df = purrr::imap_dfr(split_hpd, function(.x, .y){
    width = diff(range(.x$ObservedFrequency))
    high_sd = .x %>%
      dplyr::filter(HighSD)
    other_sd = .x %>%
      dplyr::filter(!HighSD)
    data.frame(width = width,
               highsd = median(high_sd$ObservedFrequencySD),
               othersd = median(other_sd$ObservedFrequencySD),
               allsd = median(.x$ObservedFrequencySD))
  })
  tmp_df$sample_id = hpd_res$sample_id
  tmp_df
}

plot_hpds = function(hpd_res){
  peak_data = hpd_res$peak_data
  highsd = peak_data %>%
    dplyr::filter(HighSD) %>%
    dplyr::mutate(Height = Height + 10000)
  hpd = peak_data %>%
    dplyr::filter(HPD) %>%
    dplyr::mutate(Height = Height + 100000)

  ggplot(peak_data, aes(x = ObservedMZ, xend = ObservedMZ, y = 0, yend = log10(Height))) +
    geom_segment() +
    geom_segment(data = highsd, aes(y = 6), color = "red") +
    geom_segment(data = hpd, aes(y = 7), color = "blue")
}


combine_hpd_width = function(...){
  all_width = list(...)
  purrr::map_df(all_width, ~ .x)
}
