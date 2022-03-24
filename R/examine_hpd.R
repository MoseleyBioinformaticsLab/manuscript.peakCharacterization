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


hpds_from_excel = function(scanlevel_99, scanlevel_00, msnbase, excel_files){
  sample = scanlevel_99$char_obj$zip_ms$id
  match_excel = grep(sample, excel_files, value = TRUE)

  xl_data = readxl::read_excel(match_excel, skip = 8, col_names = FALSE)
  xl_data = xl_data[, 1:2]
  names(xl_data) = c("mz", "intensity")
  message("got xl data")
  frequency_coefficients = scanlevel_99$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions$metadata$frequency_coefficients
  frequency_description = scanlevel_99$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions$metadata$frequency_fit_description
  mz_coefficients = scanlevel_99$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions$metadata$mz_coefficients
  mz_description = scanlevel_99$char_obj$zip_ms$peak_finder$peak_regions$frequency_point_regions$metadata$mz_fit_description
  org_frequency_peaks = purrr::map_df(seq_len(length(scanlevel_99$char_obj$zip_ms$peak_finder$peak_regions$peak_region_list)),
                                      function(.x){
                                        tmp = as.data.frame(scanlevel_99$char_obj$zip_ms$peak_finder$peak_regions$peak_region_list[[.x]]$peaks[, c("ObservedMZ", "ObservedFrequency", "scan")])
                                        tmp$region = .x
                                        tmp
                                      })
  xl_data$frequency = FTMS.peakCharacterization:::predict_exponentials(xl_data$mz, frequency_coefficients, frequency_description)

  sliding_metadata = scanlevel_99$char_obj$zip_ms$peak_finder$peak_regions$sliding_regions@metadata

  freq_diff = 1000
  sliding_metadata$point_spacing = round(freq_diff / 10)


  frequency_range = c(min(c(min(c(xl_data$frequency, org_frequency_peaks$ObservedFrequency)))),
                      max(c(max(c(xl_data$frequency, org_frequency_peaks$ObservedFrequency)))))

  xl_hpd = density_calculate(xl_data, sliding_metadata, frequency_range)

  peak_99 = scanlevel_99$char_obj$zip_ms$peak_finder$peak_regions$peak_data %>%
    dplyr::select(ObservedMZ, ObservedFrequency, Height, HighSD) %>%
    dplyr::mutate(source = "scanlevel_99")
  peak_00 = scanlevel_00$char_obj$zip_ms$peak_finder$peak_regions$peak_data %>%
    dplyr::select(ObservedMZ, ObservedFrequency, Height, HighSD) %>%
    dplyr::mutate(source = "scanlevel_00")
  peak_msnbase = msnbase$comb %>%
    dplyr::transmute(ObservedMZ = mz,
                     ObservedFrequency = FTMS.peakCharacterization:::predict_exponentials(ObservedMZ, frequency_coefficients, frequency_description),
                     Height = intensity,
                     HighSD = FALSE,
                     source = "msnbase")
  peak_xcal = xl_data %>%
    dplyr::transmute(ObservedMZ = mz,
                     ObservedFrequency = frequency,
                     Height = intensity,
                     HighSD = FALSE,
                     source = "xcalibur")

  split_xl = split(xl_hpd$hpd_points, xl_hpd$hpd_points$region)
  xl_ranges = purrr::map(split_xl, ~ range(.x$frequency))

  inside_peaks = purrr::imap_dfr(xl_ranges, function(in_range, range_id){
    tmp_99 = peak_99 %>%
      dplyr::filter(dplyr::between(ObservedFrequency, in_range[1], in_range[2])) %>%
      dplyr::mutate(HPDID = range_id)
    tmp_00 = peak_00 %>%
      dplyr::filter(dplyr::between(ObservedFrequency, in_range[1], in_range[2])) %>%
      dplyr::mutate(HPDID = range_id)
    tmp_msn = peak_msnbase %>%
      dplyr::filter(dplyr::between(ObservedFrequency, in_range[1], in_range[2])) %>%
      dplyr::mutate(HPDID = range_id)
    tmp_xcal = peak_xcal %>%
      dplyr::filter(dplyr::between(ObservedFrequency, in_range[1], in_range[2])) %>%
      dplyr::mutate(HPDID = range_id)
    dplyr::bind_rows(tmp_99,
                     tmp_00,
                     tmp_msn,
                     tmp_xcal)
  })
  list(peak_data = dplyr::bind_rows(peak_99, peak_00, peak_msnbase, peak_xcal),
       hpd_data = inside_peaks,
       sample_id = rename_samples(sample))
}


chisq_hpds = function(hpd_res){
  peak_data = hpd_res$hpd_data %>%
    dplyr::filter(source %in% "scanlevel_99")
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
  peak_data = hpd_res$hpd_data


  n_any = peak_data %>%
    dplyr::group_by(source, HPDID) %>%
    dplyr::summarise(n = n(),
                     n_lowsd = n - sum(HighSD))
  n_wide = n_any %>%
    tidyr::pivot_wider(id_cols = HPDID, names_from = source, values_from = n, values_fill = 0) %>%
    dplyr::arrange(dplyr::desc(scanlevel_99))
  n_wide_lowsd = n_any %>%
    tidyr::pivot_wider(id_cols = HPDID, names_from = source, values_from = n_lowsd, values_fill = 0) %>%
    dplyr::arrange(dplyr::desc(scanlevel_99)) %>%
    dplyr::transmute(HPDID = HPDID,
                     scanlevel_99_lowsd = scanlevel_99)

  all_n_wide = dplyr::left_join(n_wide, n_wide_lowsd, by = "HPDID")
  all_n_long = all_n_wide %>%
    tidyr::pivot_longer(cols = !(c(HPDID, xcalibur)),
                        names_to = "method",
                        values_to = "count")


  use_hpdid = n_wide %>%
    dplyr::slice(1) %>%
    dplyr::pull(HPDID)

  peaks_use = peak_data %>%
    dplyr::filter(HPDID %in% use_hpdid)

  peaks_type = split(peaks_use, peaks_use$source)
  peaks_type$scanlevel_99_lowsd = peaks_type$scanlevel_99 %>%
    dplyr::filter(!HighSD) %>%
    dplyr::mutate(source = "scanlevel_99_lowsd")
  plot_range = range(peaks_use$ObservedFrequency)

  plots_type = purrr::map(peaks_type, function(in_type){
    mode_height = metabolomicsUtilities::calculate_mode(in_type$Height) + 10000

    base_plot = in_type %>%
      ggplot(aes(x = ObservedFrequency,
                 xend = ObservedFrequency,
                 y = 0,
                 yend = Height)) +
      geom_segment() +
      coord_cartesian(xlim = plot_range,
                      ylim = c(NA, mode_height)) +
      labs(y = "Intensity", subtitle = in_type$source[1])
    base_plot
  })
  plots_type = plots_type[c("xcalibur", "scanlevel_00",
                            "scanlevel_99",
                            "scanlevel_99_lowsd",
                            "msnbase")]

  list(n = all_n_long,
       plots = plots_type)

}


combine_hpd_width = function(...){
  all_width = list(...)
  purrr::map_df(all_width, ~ .x)
}


get_xcalibur_peaks = function(in_file){
  xl_data = readxl::read_excel(in_file, skip = 8, col_names = FALSE)
  xl_data
}
