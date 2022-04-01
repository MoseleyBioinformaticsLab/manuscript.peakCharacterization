match_assigned_peaks = function(filtersd_assign, xcalibur_data, msnbase_data){
  # filtersd_assign = tar_read(aa_filtersd_1ecf)
  # xcalibur_data = tar_read(xcalibur_1ecf)
  # msnbase_data = tar_read(msnbase_1ecf)
  xcalibur_peaks = xcalibur_data$comb
  xcalibur_peaks$PeakID = paste0("xcalibur_", seq(1, nrow(xcalibur_peaks)))
  xcalibur_peaks$rowid = seq(1, nrow(xcalibur_peaks))
  msnbase_peaks = msnbase_data$comb
  msnbase_peaks$PeakID = paste0("msnbase_", seq(1, nrow(msnbase_peaks)))
  msnbase_peaks$rowid = seq(1, nrow(msnbase_peaks))
  assigned_peaks = purrr::map_df(filtersd_assign$emfs, ~ .x)
  theoretical_peaks = assigned_peaks %>%
    dplyr::transmute(theoreticalmz = ObservedMZ - mass_error,
                     ppm = theoreticalmz * 2e-6,
                     lowmz = theoreticalmz - ppm,
                     highmz = theoreticalmz + ppm)
  theoretical_peaks$PeakID = paste0("scanlevel_", seq(1, nrow(theoretical_peaks)))
  scanlevel_peaks = assigned_peaks %>%
    dplyr::transmute(mz = ObservedMZ,
                     intensity = Height)
  scanlevel_peaks$PeakID = theoretical_peaks$PeakID
  scanlevel_peaks$rowid = seq(1, nrow(scanlevel_peaks))
  scanlevel_peaks$mz_diff = scanlevel_peaks$mz - theoretical_peaks$theoreticalmz
  scanlevel_peaks$source = "scanlevel"

  xcalibur_matched = sub_match(theoretical_peaks, xcalibur_peaks) %>%
    dplyr::filter(grepl("^scanlevel", PeakID)) %>%
    dplyr::mutate(source = "xcalibur")
  msnbase_matched = sub_match(theoretical_peaks, msnbase_peaks) %>%
    dplyr::filter(grepl("^scanlevel", PeakID)) %>%
    dplyr::mutate(source = "msnbase")

  peak_list = list(scancentric = scanlevel_peaks$PeakID,
                   xcalibur = xcalibur_matched$PeakID,
                   msnbase = msnbase_matched$PeakID)
  peak_comb_matrix = ComplexHeatmap::make_comb_mat(peak_list)
  peak_comb_matrix

  list(comb_matrix = peak_comb_matrix,
       peaks = rbind(scanlevel_peaks,
                     xcalibur_matched,
                     msnbase_matched),
       ft_table = comb_matrix_2_table(peak_comb_matrix))

}

match_unassigned_peaks = function(filtersd_data, xcalibur_data, msnbase_data){
  # filtersd_data = tar_read(method_filtersd_97lipid)
  # xcalibur_data = tar_read(xcalibur_97lipid)
  # msnbase_data = tar_read(msnbase_97lipid)

  filtersd_peaks = filtersd_data$char_obj$zip_ms$peak_finder$peak_regions$peak_data
  xcalibur_peaks = xcalibur_data$comb
  xcalibur_peaks$PeakID = paste0("xcalibur_", seq(1, nrow(xcalibur_peaks)))
  xcalibur_peaks$rowid = seq(1, nrow(xcalibur_peaks))
  msnbase_peaks = msnbase_data$comb
  msnbase_peaks$PeakID = paste0("msnbase_", seq(1, nrow(msnbase_peaks)))
  msnbase_peaks$rowid = seq(1, nrow(msnbase_peaks))

  filtersd_tomatch = filtersd_peaks %>%
    dplyr::transmute(theoreticalmz = ObservedMZ,
                     lowmz = Start,
                     highmz = Stop,
                     PeakID = paste0("scancentric_", PeakID))

  xcalibur_peaks_matched = sub_match(filtersd_tomatch, xcalibur_peaks)
  msnbase_peaks_matched = sub_match(filtersd_tomatch, msnbase_peaks)

  peak_list = list(scancentric = filtersd_tomatch$PeakID,
                   xcalibur = xcalibur_peaks_matched$PeakID,
                   msnbase = msnbase_peaks_matched$PeakID)
  peak_comb_matrix = ComplexHeatmap::make_comb_mat(peak_list)
  peak_comb_matrix

  list(comb_matrix = peak_comb_matrix,
       ft_table = comb_matrix_2_table(peak_comb_matrix))
}

sub_match = function(filtersd_tomatch, other_peaks){
  other_peaks$mz_diff = NA
  for (irow in seq_len(nrow(filtersd_tomatch))) {
    #message(irow)
    tmp_filtersd = filtersd_tomatch[irow, ]
    other_match = dplyr::between(other_peaks$mz, tmp_filtersd$lowmz, tmp_filtersd$highmz)
    if (sum(other_match) == 1) {
      other_peaks[other_match, "PeakID"] = tmp_filtersd$PeakID
      other_peaks[other_match, "mz_diff"] = other_peaks[other_match, "mz"] - tmp_filtersd$theoreticalmz
    } else if (sum(other_match) > 1) {
      tmp_other = other_peaks[other_match, ]
      min_loc = which.min(abs(tmp_other$mz - tmp_filtersd$theoreticalmz))
      other_peaks[tmp_other$rowid[min_loc], "PeakID"] = tmp_filtersd$PeakID
      other_peaks[tmp_other$rowid[min_loc], "mz_diff"] = other_peaks[tmp_other$rowid[min_loc], "mz"] - tmp_filtersd$theoreticalmz
    }
  }
  other_peaks
}

comb_matrix_2_table = function(comb_matrix){
  df_comb = as.data.frame(comb_matrix)
  org_vars = names(df_comb)
  set_sizes = ComplexHeatmap::set_size(comb_matrix)
  comb_sizes = ComplexHeatmap::comb_size(comb_matrix)

  df_comb = cbind(df_comb, set_sizes)
  df_comb = rbind(df_comb, c(comb_sizes, NA))
  df_samples = data.frame(method = rownames(df_comb))
  df_comb = cbind(df_samples, df_comb)
  rownames(df_comb) = NULL

  limit_col = ncol(comb_matrix)
  limit_row = nrow(comb_matrix)

  df_comb2 = purrr::imap_dfc(df_comb, function(in_col, col_id){
    out_col = as.character(in_col)
    row_index = seq(1, length(in_col))
    if (col_id %in% "method") {
      out_col[row_index > limit_row] = ""
    }
    if (col_id %in% "set_sizes") {
      out_col[row_index > limit_row] = ""
    }
    if (col_id %in% org_vars) {
      out_col[(row_index <= limit_row) & (out_col %in% "1")] = "x"
      out_col[(row_index <= limit_row) & (out_col %in% "0")] = ""
    }
    out_col
  })

  df_comb2$method[nrow(df_comb2)] = "comb_sizes"

  out_ft = flextable(df_comb2)

  new_labels = names(df_comb2)
  new_labels[new_labels %in% org_vars] = ""
  names(new_labels) = names(df_comb2)
  out_ft = set_header_labels(out_ft, values = new_labels)
  out_ft = align(out_ft, i = nrow(df_comb2), j = seq(2, ncol(df_comb2)), align = "right", part = "body")
  out_ft = align(out_ft, i = seq(1, (nrow(df_comb2) - 1)), j = seq(2, (ncol(df_comb2) - 1)), align = "center", part = "body")
  # out_ft = hline(out_ft, i = nrow(df_comb2) - 1, part = "body")
  # out_ft = vline(out_ft, j = ncol(df_comb2) - 1, part = "body")
  out_ft
}
