# create region counts
create_regions_count_histogram <- function(in_char, file_loc){
  peak_regions <- in_char$peak_finder$peak_regions

  counted_regions <- FTMS.peakCharacterization:::count_overlaps(peak_regions$sliding_regions, peak_regions$mz_point_regions)
  count_data <- as.data.frame(counted_regions@elementMetadata)

  p <- ggplot(count_data, aes(x = nonzero_counts)) + geom_histogram(bins = 100) +
    geom_vline(xintercept = quantile(count_data$nonzero_counts, .99), color = "red") + scale_y_log10(limits = c(1, NA), expand = c(0,0)) + labs(y = "Log10(count)", x = "# non-zero counts")
  base_file <- tools::file_path_sans_ext(file_loc)
  svg_file <- paste0(base_file, ".svg")
  Cairo::CairoSVG(file = svg_file, width = 8, height = 5)
  print(p)
  dev.off()
  sys_call <- paste0("inkscape -z -d 300 -e ", file_loc, " ", svg_file)
  system(sys_call)
}

# find and plot a multiple peak region
#
# We will cheat, and do reduction again, and then take the fully split out,
# and find those ones where there were originally multiples in the region
create_multiple_peak_figure <- function(characterized_sample, basic_sample, file_loc){
  characterized_regions <- characterized_sample$peak_regions$clone(deep = TRUE)
  sample_regions <- basic_sample$peak_finder$clone(deep = TRUE)

  sample_regions$reduce_sliding_regions()

  characterized_peaks <- characterized_regions$peak_data %>% mutate(mz = ObservedMZ)
  characterized_peaks <- mz_points_to_regions(characterized_peaks)

  sample_init_regions <- sample_regions$peak_regions$peak_regions

  num_overlaps <- IRanges::countOverlaps(sample_init_regions, characterized_peaks)

  use_region <- sample_init_regions[which.max(num_overlaps)]

  mz_point_regions <- IRanges::subsetByOverlaps(characterized_regions$mz_point_regions, use_region)
  point_data <- as.data.frame(mz_point_regions@elementMetadata) %>% mutate(scan = as.factor(scan))

  point_plot <- ggplot(point_data, aes(x = mz, y = log10(intensity), color = scan)) +
    geom_point() + geom_line() + theme(legend.position = "none") + labs(x = "M/Z", y = "Log10(Intensity)")

  tiled_regions <- IRanges::subsetByOverlaps(sample_regions$peak_regions$tiled_regions, use_region)

  mz_point_regions@elementMetadata$log_int <- log(mz_point_regions@elementMetadata$intensity + 1e-8)
  mz_point_regions <- split(mz_point_regions, mz_point_regions@elementMetadata$scan)

  reduced_peaks <- purrr::map_df(names(mz_point_regions), function(in_scan){
    FTMS.peakCharacterization:::get_reduced_peaks(mz_point_regions[[in_scan]], peak_method = "lm_weighted", min_points = 4)
  })

  reduced_peaks <- reduced_peaks[!is.na(reduced_peaks$ObservedMZ), ] %>% mutate(scan = as.factor(scan),
                                                                                mz = ObservedMZ)

  reduced_mz_points <- mz_points_to_regions(reduced_peaks, mz_point_regions[[1]]@metadata$point_multiplier)
  tiled_regions@elementMetadata$peak_count <- IRanges::countOverlaps(tiled_regions, reduced_mz_points)
  tiled_regions_data <- as.data.frame(tiled_regions@elementMetadata)
  tiled_regions_data <- mutate(tiled_regions_data, mz = (mz_end + mz_start) / 2,
                               width = mean(mz_end - mz_start), scan = as.factor(1))
  count_peak_intensity_ratio <- max(log10(reduced_peaks$Height)) / max(tiled_regions_data$peak_count)

  tiled_regions_data <- mutate(tiled_regions_data, scaled_counts = count_peak_intensity_ratio * peak_count)

  tiled_segments <- IRanges::reduce(tiled_regions[tiled_regions_data$peak_count > 0])
  point_factor <- reduced_mz_points@metadata$point_multiplier
  tiled_segments_mz <- data.frame(mz_start = IRanges::start(tiled_segments) / point_factor,
                                  mz_end = IRanges::end(tiled_segments) / point_factor,
                                  count = 2)

  reduced_peak_plot <- ggplot(reduced_peaks, aes(x = ObservedMZ, y = log10(Height), color = scan)) + geom_point(color = "gray") + theme(legend.position = "none") + geom_point(data = tiled_regions_data, aes(x = mz, y = scaled_counts), color = "black") + scale_y_continuous(sec.axis = sec_axis(~./count_peak_intensity_ratio, name = "# of Scan Peaks")) + geom_segment(data = tiled_segments_mz, aes(x = mz_start, xend = mz_end, y = count, yend = count), color = "blue", size = 2)


  out_plot <- point_plot + reduced_peak_plot + plot_annotation(tag_levels = "A")

  base_file <- tools::file_path_sans_ext(file_loc)
  svg_file <- paste0(base_file, ".svg")
  Cairo::CairoSVG(file = svg_file, width = 16, height = 5)
  print(out_plot)
  dev.off()
  sys_call <- paste0("inkscape -z -d 300 -e ", file_loc, " ", svg_file)
  system(sys_call)
}
