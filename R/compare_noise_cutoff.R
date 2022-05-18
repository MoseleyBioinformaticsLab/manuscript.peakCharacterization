#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param in_char
#' @return
#' @author rmflight
#' @export
compare_noise_cutoffs <- function(in_char) {

  use_percentiles = c(seq(0.1, 0.9, by = 0.1),
                      seq(0.91, 0.99, by = 0.01))

  peak_finder = in_char$zip_ms$peak_finder$clone(deep = TRUE)

  n_region = purrr::map_df(use_percentiles, function(in_percentile){
    out_regions = find_signal_adj(peak_finder$peak_regions$sliding_regions,
                                  peak_finder$peak_regions$frequency_point_regions$frequency,
                                  peak_finder$quantile_multiplier,
                                  peak_finder$n_point_region,
                                  in_percentile)
    data.frame(percentile = in_percentile,
               n_region = length(out_regions),
               sample = in_char$zip_ms$id
               )
  })

  zero_regions = reduce_removing_zero(peak_finder$peak_regions$sliding_regions,
                                      peak_finder$peak_regions$frequency_point_regions$frequency,
                                      peak_finder$quantile_multiplier,
                                      peak_finder$n_point_region)

  n_region = rbind(n_region,
                   data.frame(percentile = 0,
                              n_region = length(zero_regions),
                              sample = in_char$zip_ms$id))

  n_region$sample = rename_samples(n_region$sample)
  n_region
}

plot_example_sliding_counts = function(peak_data){
  #source("./packages.R")
  #peak_data = tar_read(data_97lipid)
  sliding_regions = peak_data$zip_ms$peak_finder$peak_regions$sliding_regions
  point_region_list = peak_data$zip_ms$peak_finder$peak_regions$frequency_point_regions$frequency
  region_multiplier = peak_data$zip_ms$peak_finder$peak_regions$frequency_multiplier
  multiplier = peak_data$zip_ms$peak_finder$quantile_multiplier
  n_point_region = peak_data$zip_ms$peak_finder$n_point_region

  nz_counts = FTMS.peakCharacterization:::count_overlaps(sliding_regions, point_region_list[[1]])
  n_region = seq(2, length(point_region_list))
  for (iregion in n_region) {
    nz_counts_iregion = FTMS.peakCharacterization:::count_overlaps(sliding_regions, point_region_list[[iregion]])
    nz_counts = nz_counts + nz_counts_iregion
  }

  mcols(sliding_regions) = DataFrame(nz_counts = nz_counts, frequency = start(sliding_regions) / region_multiplier)

  chunk_indices = seq(1, length(nz_counts), by = n_point_region)

  chunk_perc = purrr::map_dbl(chunk_indices, function(in_index){
    use_counts = nz_counts[seq(in_index, min(in_index + (n_point_region - 1), length(nz_counts)))]
    if (max(use_counts) > 0) {
      return(stats::quantile(use_counts, 0.99))
    } else {
      return(0)
    }
  })

  chunk_cutoff = ceiling(median(chunk_perc) * multiplier)

  freq_800 = FTMS.peakCharacterization:::predict_exponentials(800, peak_data$peak_finder$peak_regions$frequency_point_regions$metadata$frequency_coefficients, peak_data$peak_finder$peak_regions$frequency_point_regions$metadata$frequency_fit_description)[1]
  iranges_800 = IRanges::IRanges(start = (freq_800 - 500) * 400, end = (freq_800 + 500) * 400)

  point_list_nz = purrr::map(point_region_list, function(.x){
    tmp_df = as.data.frame(mcols(.x))
    nz = tmp_df$intensity > 0
    .x[nz]
  })

  point_list_800 = as.data.frame(mcols(purrr::map(point_list_nz, ~ IRanges::subsetByOverlaps(.x, iranges_800)) %>% IRanges::IRangesList() %>% unlist()))

  sliding_800 = IRanges::subsetByOverlaps(sliding_regions, iranges_800)
  sliding_800_data = as.data.frame(mcols(sliding_800))
  sliding_800_data$nz_plot = sliding_800_data$nz_counts
  sliding_800_data$nz_plot[sliding_800_data$nz_counts == 0] = 1
  sliding_800_data_99 = stats::quantile(sliding_800_data$nz_counts, 0.99)

  xlim = c(start(iranges_800), end(iranges_800)) / 400

  keep_sliding_800 = sliding_800_data$nz_counts > chunk_cutoff

  sliding_800_data_keep = sliding_800_data[keep_sliding_800, ]

  y_title_size = 7.5


  intensity_plot = ggplot(point_list_800, aes(x = frequency, y = log10(intensity))) +
    geom_point(alpha = 0.5) +
    coord_cartesian(xlim = xlim, ylim = c(2, NA)) +
    labs(y = "Log10(Intensity)", x = "") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    theme(axis.title.y = element_text(size = y_title_size))


  sliding_plot = ggplot(sliding_800_data, aes(x = frequency, y = nz_counts)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = chunk_cutoff, color = "red") +
    geom_hline(yintercept = sliding_800_data_99, color = "blue") +
    coord_cartesian(xlim = xlim) +
    labs(y = "Non-Zero Counts", x = "") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    theme(axis.title.y = element_text(size = y_title_size))


  sliding_plot_zoom = sliding_plot + ylim(0, 10)

  keep_sliding_plot = ggplot(sliding_800_data_keep, aes(x = frequency, y = nz_counts)) +
    geom_point() +
    coord_cartesian(xlim = xlim) +
    labs(y = "Non-Zero Counts", x = "") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    theme(axis.title.y = element_text(size = y_title_size))




  # this looks exactly like I expect! Awesome
  #intensity_plot / sliding_plot / sliding_plot_zoom / keep_sliding_plot

  sliding_800_keep = reduce(sliding_800[keep_sliding_800])

  point_list_keep_800 = as.data.frame(mcols(purrr::map(point_list_nz, ~ IRanges::subsetByOverlaps(.x, sliding_800_keep)) %>% IRanges::IRangesList() %>% unlist()))

  intensity_keep_plot = ggplot(point_list_keep_800, aes(x = frequency, y = log10(intensity))) +
    geom_point(alpha = 0.5) +
    coord_cartesian(xlim = xlim, ylim = c(2, NA)) +
    labs(y = "Log10(Intensity)", x = "Frequency") +
    theme(axis.title.y = element_text(size = y_title_size))


  out_plot = (intensity_plot / sliding_plot / sliding_plot_zoom / keep_sliding_plot / intensity_keep_plot)

  out_plot
}

n_final_peak = function(char_list){
  char_obj = char_list$char_obj
  n_peak = nrow(char_obj$zip_ms$peak_finder$peak_regions$peak_data)
  data.frame(processed = char_list$processed,
             n_peak = n_peak,
             sample = rename_samples(char_obj$zip_ms$id))
}

plot_sliding_window_density = function(peak_data){
  # peak_data = tar_read(data_97lipid)
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

  single_chunk_region = sliding_regions[1:n_point_region]
  single_region = sliding_regions[1]

  single_chunk = data.frame(use_counts = nz_counts[
    seq(chunk_indices[1], (chunk_indices[2] - 1))])
  single_99 = stats::quantile(single_chunk$use_counts, 0.99)

  all_99 = data.frame(perc_99 = chunk_perc,
                      index = seq(1, length(chunk_perc)))
  all_cut = ceiling(1.5 * median(all_99$perc_99))
  single_plot = ggplot(single_chunk, aes(x = use_counts)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = single_99, color = "blue") +
    scale_y_log10(expand = c(0, 0), limits = c(1, NA)) +
    labs(x = "Number of Non-Zero Points",
         y = "Log10(Counts)")
  all_plot = ggplot(all_99, aes(x = perc_99)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = all_cut, color = "red") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = "# Tiled Regions",
         x = "99th Percentiles of Tiled Regions")
  all_plot
}


create_noise_plot = function(noise_combine){
  ecf_plot = noise_combine %>%
    dplyr::filter(grepl("ecf", sample)) %>%
      ggplot(aes(x = percentile, y = n_region, color = sample)) +
    geom_point() +
    geom_line() +
    theme(legend.position = c(0.1, 0.5)) +
    labs(y = "# of Initial Peak Regions", x = "Percentile Cutoff")

  lipid_plot = noise_combine %>%
    dplyr::filter(grepl("lipid", sample)) %>%
    ggplot(aes(x = percentile, y = n_region, color = sample)) +
    geom_point() +
    geom_line() +
    theme(legend.position = c(0.1, 0.5)) +
    labs(y = "# of Initial Peak Regions", x = "Percentile Cutoff")

  (ecf_plot | lipid_plot)
}

calculate_m_n_ratio = function(noise_combine){
  noise_by_sample = split(noise_combine, noise_combine$sample)
  mn_sample = purrr::map(noise_by_sample, function(in_noise){
    count_0 = in_noise[ in_noise$percentile == 0, "n_region"]
    count_99 = in_noise[ in_noise$percentile == 0.99, "n_region"]
    ratio_val = english::as.english(round(count_0 / count_99))

    list(count_0 = count_0,
         count_99 = count_99,
         ratio = ratio_val)

  })
  mn_sample
}
