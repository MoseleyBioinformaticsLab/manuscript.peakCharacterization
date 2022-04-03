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
  if (grepl("xcalibur", other_file)) {
    other_data = purrr::map(other_data, function(.x){
      .x %>%
        dplyr::transmute(intensity = Intensity,
                         mz = TargetMass,
                         sample_id = sample)
    })
  }
  if (grepl("msnbase", other_file)) {
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

  na_vector = rep(NA, nrow(scancentric_location))
  names(na_vector) = rownames(scancentric_location)

  split_location = purrr::map(colnames(scancentric_location), function(isample){
    scancentric_location[, isample]
  })
  names(split_location) = colnames(scancentric_location)

  match_loc_intensity = furrr::future_imap(split_location, function(in_location, isample){
    out_location = out_intensity = na_vector
    location_ppm = in_location * 2e-6
    location_low = in_location - location_ppm
    location_high = in_location + location_ppm

    tmp_df = other_data[[isample]]
    for (imf in names(in_location)) {
      if (!is.na(in_location[imf])) {
        match_locs = tmp_df %>%
          dplyr::filter(dplyr::between(mz, location_low[imf], location_high[imf])) %>%
          dplyr::mutate(mz_diff = abs(mz - in_location[imf])) %>%
          dplyr::arrange(dplyr::desc(mz_diff))

        if (nrow(match_locs) > 0) {
          out_location[imf] = match_locs$mz[1]
          out_intensity[imf] = match_locs$intensity[1]
        }
      }
    }
    list(location = out_location,
         intensity = out_intensity)
  })

  other_locations = as.matrix(purrr::map_dfc(match_loc_intensity, ~ .x$location))
  rownames(other_locations) = rownames(scancentric_location)
  other_intensity = as.matrix(purrr::map_dfc(match_loc_intensity, ~ .x$intensity))
  rownames(other_intensity) = rownames(scancentric_location)

  other_imfs = scancentric_imfs
  other_imfs$intensity = other_intensity
  other_imfs$location = other_location
  other_imfs
}
