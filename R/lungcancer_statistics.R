create_sample_info_df = function(sample_list_file){
  sample_info = read.table(sample_list_file, header = TRUE, sep = "\t")
  sample_info
}

extract_scancentric_imfs = function(emf_file){
  emf_data = readRDS(emf_file)
  imf_height = smirfeTools::extract_imf_emf_data(emf_data, by = "IMF", intensity = Height)
  imf_corrected_height = smirfeTools::extract_imf_emf_data(emf_data, by = "IMF", intensity = CorrectedLog10Height)

  raw_intensity = imf_height$intensity
  cor_intensity = 10^imf_corrected_height$intensity

  imf_height$corrected_intensity = cor_intensity
  imf_height
}

match_imfs = function(other_file, scancentric_imfs){
  other_data = readRDS(other_file)
  if (grepl("xcalibur")) {
    other_data = purrr::map(other_data, function(.x){
      dplyr::transmute(intensity = Intensity,
                       mz = TargetMass,
                       sample_id = sample)
    })
  }
  if (grepl("msnbase")) {
    other_data = purrr::map(other_data, function(.x){
      tmp_df = .x$comb %>%
        dplyr::transmute(intensity = intensity,
                         mz = mz,
                         sample_id = .x$sample_id)
    })
  }
  other_samples = purrr::map_chr(other_data, ~ .x$sample_id[1])
  names(other_data) = other_samples
  scancentric_location = scancentric_imfs$location
  keep_samples = intersect(other_samples, colnames(scancentric_location))
  keep_other = other_samples %in% keep_samples
  other_data = other_data[keep_other]
  other_data = other_data[colnames(scancentric_location)]

  other_location = other_intensity = matrix(NA, nrow = nrow(scancentric_location),
                                            ncol = ncol(scancentric_location))
  rownames(other_location) = rownames(other_intensity) = rownames(scancentric_location)
  colnames(other_location) = colnames(other_intensity) = colnames(scancentric_location)
  location_ppm = scancentric_location * 2e-6
  location_low = scancentric_location - location_ppm
  location_high = scancentric_location + location_ppm

  for (icol in colnames(scancentric_location)) {
    tmp_df = other_data[[icol]] %>%
      dplyr::mutate(irow = seq(1, nrow(.)))

    for (irow in rownames(scancentric_location)) {
      if (!is.na(scancentric_location[irow, icol])) {
        match_locs = dplyr::between(tmp_df$mz, location_low[irow, icol], location_high[irow, icol])
      }
    }

  }
}
