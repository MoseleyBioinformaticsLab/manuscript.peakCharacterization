msnbase_only = function(mzml_file){
  mzml_prof = MSnbase::readMSData(mzml_file, msLevel. = 1, centroided. = FALSE)
  mzml_info = get_ms_info(mzml_prof)

  all_scans_cent = mzml_prof %>%
    MSnbase::pickPeaks()
  all_scans_cent_mz = mz(all_scans_cent)
  all_scans_cent_intensity = intensity(all_scans_cent)

  all_scans_data = purrr::map_df(seq_len(length(all_scans_cent_mz)), function(in_scan){
    data.frame(mz = all_scans_cent_mz[[in_scan]],
               intensity = all_scans_cent_intensity[[in_scan]],
               scan = in_scan)
  })

  comb_prof = MSnbase::combineSpectra(mzml_prof, method = meanMzInts, mzd = 0, ppm = 1)

  comb_cent = comb_prof %>%
    MSnbase::pickPeaks()
  comb_cent_data = data.frame(mz = mz(comb_cent)[[1]],
                              intensity = intensity(comb_cent)[[1]])

  list(scanlevel = all_scans_data, comb = comb_cent_data)
}
