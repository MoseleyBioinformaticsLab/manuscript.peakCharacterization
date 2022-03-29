get_assignments = function(assignments_file){
  #assignments_file = tar_read(assign_filtersd_1ecf)
  assign_data = smirfeTools::read_smirfe_assignment(assignments_file$file, assigned_only = FALSE)
  #score_filter_assignments(assign_data, filter_conditions = e_value <= 0.5)
}

get_expected_formulas = function(){
  formulas = read.table("data/data_input/expected_formulas.txt",
                        header = TRUE, sep = ":")
  formulas
}

check_nap_intensity = function(nap_df){
  nap_int_ratios = calculate_log_ratio_differences(nap_df)
}

find_aa_assignments = function(assign_data, aa_formulas, e_cutoff = NULL){
  if (is.null(e_cutoff)) {
    scored_assignments = smirfeTools::score_filter_assignments(assign_data)
  } else {
    scored_assignments = smirfeTools::score_filter_assignments(assign_data, filter_conditions = e_value <= e_cutoff)
  }

  within_emfs = smirfeTools:::get_sample_emfs(scored_assignments$assignments, assign_data$sample, evalue_cutoff = 0.98, use_corroborating = TRUE)

  possible_aa = purrr::map_lgl(within_emfs, function(in_emf){
    if (any(aa_formulas$complete_EMF %in% in_emf$complete_EMF)) {
      TRUE
    } else {
      FALSE
    }
  })
  within_aa = within_emfs[possible_aa]
  all_aa = purrr::map_df(within_aa, ~ .x)
  all_aa = dplyr::left_join(all_aa, dplyr::select(aa_formulas, -isotopologue_EMF), by = "complete_EMF")
  all_aa = dplyr::left_join(dplyr::select(all_aa, -PeakID, -Sample, -ObservedMZ), scored_assignments$data, by = "Sample_Peak")

  split_aa = split(all_aa, all_aa$AA)
  list(aa = split_aa,
       aa_formulas = aa_formulas,
       scored = scored_assignments)
}

get_number_sample_peak = function(sample_str){
  # sample_str = "161212_unlabeledAAs_1_ECF_1805"
  split_str = strsplit(sample_str, "_", fixed = TRUE)
  num_id = purrr::map_chr(split_str, function(.x){
    .x[length(.x)]
  })

  num_id
}

#' log ratio differences
#'
#' @param nap_df
#'
#' @details nap_df should be a data.frame containing:
#'   NAP: log10 transformed NAP values
#'   intensity: log10 transformed intensity values
#'   Sample_Peak: the peak identifier
calculate_log_ratio_differences = function(nap_df){
  n_comp = nrow(nap_df) * (nrow(nap_df) - 1) / 2
  comp_df = data.frame(p1 = rep(NA, n_comp),
                       p2 = rep(NA, n_comp),
                       nap_ratio = rep(NA, n_comp),
                       intensity_ratio = rep(NA, n_comp))
  comp_count = 1
  for (i in seq_len(nrow(nap_df) - 1)) {

    for (j in seq(i + 1, nrow(nap_df))) {
      #message(paste0(i, "_", j))
      comp_df[comp_count, "nap_ratio"] = (nap_df$NAP[i]) - (nap_df$NAP[j])
      comp_df[comp_count, "intensity_ratio"] = (nap_df$intensity[i]) - (nap_df$intensity[j])
      comp_df[comp_count, "p1"] = nap_df$Sample_Peak[i]
      comp_df[comp_count, "p2"] = nap_df$Sample_Peak[j]
      comp_count = comp_count + 1
    }

  }
  comp_df = comp_df %>%
    dplyr::mutate(nap_intensity_diff = abs(nap_ratio - intensity_ratio))
  p1 = get_number_sample_peak(comp_df$p1)
  p2 = get_number_sample_peak(comp_df$p2)

  out_id = rep(NA, n_comp)
  for (i in seq_len(n_comp)) {
    out_id[i] = paste(sort(c(p1[i], p2[i])), collapse = ":")
  }
  comp_df$p1p2 = out_id
  comp_df
}

process_imf = function(imf_string){
  element_pieces = stringr::str_extract_all(imf_string, '(\\d+[A-Z][a-z]?)(\\d+)')[[1]]
  isotopes = unlist(stringr::str_extract_all(element_pieces, '\\d+[A-Z][a-z]?'))
  counts = stringr::str_remove_all(element_pieces, '\\d+[A-Z][a-z]?')
  data.frame(isotopes = isotopes, count = counts)
}


check_serine = function(filtersd){
  # filtersd = tar_read(aa_filtersd_1ecf)
  #
  use_aa = filtersd$aa$Serine
  aa_data = dplyr::left_join(dplyr::select(use_aa, -PeakID, -Sample, -ObservedMZ), filtersd$scored$data, by = "Sample_Peak")

  tmp_data = aa_data %>%
    dplyr::select(Sample_Peak, NAP, Height, e_value, isotopologue_IMF, lbl.type, lbl.count, ObservedMZ, CorrectedLog10Height)
  tmp_data$AltNAP = tmp_data$NAP
  tmp_data$AltNAP[1] = 1

  nap_df = tmp_data %>%
    dplyr::mutate(intensity = 10^CorrectedLog10Height,
                  NAP = NAP)

  nap_ratio = calculate_log_ratio_differences(nap_df)

  nap_org = nap_df
  nap_org$NAP = tmp_data$NAP
  nap_org_ratio = calculate_log_ratio_differences(nap_org)
  write.table(tmp_data, file = "serine_example.csv", sep = ",",
              col.names = TRUE, row.names = FALSE)
  write.table(nap_ratio, file = "")

}

calculate_comparisons_matrix = function(in_matrix){
  n_sample = ncol(in_matrix)
  peaks_comp = colnames(in_matrix)
  n_comp = n_sample * (n_sample - 1) / 2
  comp_matrix = matrix(NA, nrow = nrow(in_matrix), ncol = n_comp)
  comp_df = data.frame(p1 = rep(NA, n_comp),
                       p2 = rep(NA, n_comp))
  i_comp = 1
  for (i in seq_len(n_sample - 1)) {
    for (j in seq(i + 1, n_sample)) {
      comp_matrix[, i_comp] = in_matrix[, i] - in_matrix[, j]
      comp_df$p1[i_comp] = peaks_comp[i]
      comp_df$p2[i_comp] = peaks_comp[j]

      i_comp = i_comp + 1
    }
  }
  list(matrix = comp_matrix, peaks = comp_df)
}

aa_height_nap_scanlevel = function(nap_df, scan_level){
  # aa_data = tar_read(aa_filtersd_1ecf)

  intensity_use = scan_level[nap_df$Sample_Peak]
  intensity_matrix = as.matrix(dplyr::bind_rows(intensity_use))

  nap_matrix = matrix(nap_df$NAP,
                      nrow = nrow(intensity_matrix),
                      ncol = ncol(intensity_matrix),
                      byrow = TRUE)
  colnames(nap_matrix) = colnames(intensity_matrix)

  nap_comparisons = calculate_comparisons_matrix(nap_matrix)
  intensity_comparisons = calculate_comparisons_matrix(intensity_matrix)

  nap_intensity = abs(nap_comparisons$matrix - intensity_comparisons$matrix)
  rownames(nap_intensity) = seq(1, nrow(nap_intensity))

  list(nap = nap_comparisons$matrix,
       intensity = intensity_comparisons$matrix,
       nap_intensity = nap_intensity,
       peaks = nap_comparisons$peaks)
}

match_peaks = function(theoretical_peaks, other_peaks){
  # theoretical_peaks = complete_theory
  # other_peaks = xcalibur_peaks
  matched = purrr::map_df(seq(1, nrow(theoretical_peaks)), function(in_row){
    possible_match = other_peaks %>%
      dplyr::filter(dplyr::between(mz, theoretical_peaks[in_row, "lowermz"],
                                   theoretical_peaks[in_row, "highmz"]))
    if (nrow(possible_match) == 0) {
      return(NULL)
    } else if (nrow(possible_match) == 1) {
      possible_match = possible_match[which.min(abs(possible_match$mz - theoretical_peaks[in_row, "theoreticalmz"])), ]
    }
    possible_match %>%
      dplyr::transmute(ObservedMZ = mz,
                       NAP = theoretical_peaks$NAP[in_row],
                       intensity = intensity,
                       Their_Peak = Sample_Peak,
                       Sample_Peak = theoretical_peaks$Sample_Peak[in_row],
                       isotopologue_IMF = theoretical_peaks$isotopologue_IMF[in_row])
  })
  matched
}

aa_height_nap_all = function(filtersd, xcalibur, msnbase){
  # filtersd = tar_read(aa_filtersd_1ecf)
  # xcalibur = tar_read(xcalibur_1ecf)
  # msnbase = tar_read(msnbase_1ecf)
  scan_level_data = filtersd$scored$scan_level
  xcalibur_peaks = xcalibur$comb
  xcalibur_peaks$Sample_Peak = paste0("xcalibur_", seq_len(nrow(xcalibur_peaks)))

  msnbase_peaks = msnbase$comb
  msnbase_peaks$Sample_Peak = paste0("msnbase_", seq_len(nrow(msnbase_peaks)))

  aa_info = purrr::imap(filtersd$aa, function(use_aa, aa_id){
    # use_aa = filtersd$aa[["Alanine"]]
    # aa_id = "Alanine"
    #message(aa_id)

    use_aa = use_aa %>%
      dplyr::mutate(complete_label = paste0(complete_EMF, ".", lbl.type, ".", lbl.count))

    aa_complete = split(use_aa, use_aa$complete_label)
    aa_hasenough = purrr::map_lgl(aa_complete, ~ (nrow(.x) > 1))
    aa_complete = aa_complete[aa_hasenough]


    aa_ratios = purrr::imap(aa_complete, function(in_complete, split_id){
      #message(split_id)
      use_intensities = c("CorrectedLog10Height",
                      "Log10Height")

      char_level = purrr::map(use_intensities, function(in_intensity){
        tmp_df = in_complete %>%
          dplyr::transmute(NAP = log10(NAP),
                          intensity = .data[[in_intensity]],
                          Sample_Peak = Sample_Peak)

        peak_ratios = calculate_log_ratio_differences(tmp_df)
        peak_ratios$source = in_intensity
        peak_ratios
      })
      names(char_level) = use_intensities

      scan_level = purrr::map(use_intensities, function(in_intensity){
        tmp_df = in_complete %>%
          dplyr::transmute(NAP = log10(NAP),
                           Sample_Peak = Sample_Peak)
        scan_ratios = aa_height_nap_scanlevel(tmp_df, scan_level_data[[in_intensity]])
        scan_ratios$source = in_intensity
        scan_ratios
      })
      names(scan_level) = use_intensities

      complete_theory = in_complete %>%
        dplyr::transmute(NAP = log10(NAP),
                         theoreticalmz = ObservedMZ - mass_error,
                         Sample_Peak = Sample_Peak,
                         isotopologue_IMF = isotopologue_IMF,
                         lowermz = theoreticalmz - (theoreticalmz * 2e-6),
                         highmz = theoreticalmz + (theoreticalmz * 2e-6))



      xcal_matched = match_peaks(complete_theory, xcalibur_peaks)


      if (nrow(xcal_matched) > 1) {
        xcal_matched$intensity = log10(xcal_matched$intensity)
        xcal_ratios = calculate_log_ratio_differences(xcal_matched[, c("Sample_Peak", "NAP", "intensity")])
        xcal_ratios$source = "xcalibur"
      } else {
        xcal_ratios = NULL
      }

      msnbase_matched = match_peaks(complete_theory, msnbase_peaks)
      if (nrow(msnbase_matched) > 1) {
        msnbase_matched$intensity = log10(msnbase_matched$intensity)
        msnbase_ratios = calculate_log_ratio_differences(msnbase_matched[, c("Sample_Peak", "NAP", "intensity")])
        msnbase_ratios$source = "msnbase"

      } else {
        msnbase_ratios = NULL
      }

      list(characterized = list(ratios = char_level,
                                peaks = in_complete),
           scanlevel = list(ratios = scan_level,
                            peaks = in_complete),
           xcalibur = list(ratios = xcal_ratios,
                           peaks = xcal_matched),
           msnbase = list(ratios = msnbase_ratios,
                          peaks = msnbase_matched))

    })
  })
  aa_info
}

scan_level_10th = function(scan_level_data, which = c("complete")){
  # aa_ratios = tar_read(nap_height_1ecf)
  # scan_level_data = aa_ratios$Alanine$C8H15N1Na1O4.15N.0$scanlevel$ratios$Log10Height
  ratio_diffs = scan_level_data$nap_intensity
  n_na = apply(ratio_diffs, 1, function(.x){
    sum(is.na(.x))
  })
  median_diffs = rowMedians(ratio_diffs, na.rm = TRUE)

  ratio_nona = ratio_diffs[n_na == 0, , drop = FALSE]
  median_nona = median_diffs[n_na == 0]

  median_order = order(median_nona, decreasing = FALSE)
  ratio_nona = ratio_nona[median_order, , drop = FALSE]
  median_nona = median_nona[median_order]

  n_scan = nrow(ratio_diffs)
  scan_10th = round(0.1 * n_scan)

  if (nrow(ratio_nona) >= (2 * scan_10th)) {
    return(ratio_nona[seq(1, scan_10th), , drop = FALSE])
  } else {
    return(NULL)
  }


}

diff_vector_2_df = function(diffs, source = "none"){
  out_df = data.frame(ratio = seq_len(length(diffs)),
             diff = diffs) %>%
    dplyr::mutate(source = source)
  out_df
}

diff_matrix_2_df = function(diffs, source = "none"){
  ids = seq_len(ncol(diffs))
  scans = rownames(diffs)
  diff_df = as.data.frame(diffs)
  names(diff_df) = ids
  diff_df$scan = scans
  diff_df_out = tidyr::pivot_longer(diff_df, cols = !scan, names_to = "ratio",
                      values_to = "diff") %>%
    dplyr::mutate(source = source)

  diff_df_out
}

count_missing_scans = function(scan_level_data){
  # scan_level_data = aa_formula$scanlevel$ratios$CorrectedLog10Height
  scan_level_ratios = scan_level_data$nap_intensity
  n_missing = apply(scan_level_ratios,
                    2,
                    function(.x){
    sum(is.na(.x))
                    })
  n_missing
}

compare_peak_ratios = function(aa_ratios){
  # aa_ratios = tar_read(nap_height_1ecf)

  purrr::imap(aa_ratios, function(in_aa, aa_id){
    # in_aa = aa_ratios[[1]]
    message(aa_id)
    purrr::imap(in_aa, function(aa_formula, formula_id){
      # aa_formula = aa_ratios$Glycine$C7H13N1Na1O4.15N.0
      # aa_formula = in_aa[[1]]
      message(formula_id)
      n_missing = count_missing_scans(aa_formula$scanlevel$ratios$CorrectedLog10Height) %>%
        diff_vector_2_df() %>%
        dplyr::transmute(ratio = ratio,
                         n_missing = diff)
      corrected_scan_10 = scan_level_10th(aa_formula$scanlevel$ratios$CorrectedLog10Height)

      if (is.null(corrected_scan_10)) {
        return(NULL)
      }
      corrected_scan10_max = apply(corrected_scan_10, 2, max) %>%
        diff_vector_2_df(source = "scan_max_corrected")

      corrected_char = aa_formula$characterized$ratios$CorrectedLog10Height$nap_intensity_diff %>%
        diff_vector_2_df(source = "char_corrected")

      raw_scan_10 = scan_level_10th(aa_formula$scanlevel$ratios$Log10Height)
      raw_scan10_max = apply(raw_scan_10, 2, max) %>%
        diff_vector_2_df(source = "scan_max_raw")
      raw_char = aa_formula$characterized$ratios$Log10Height$nap_intensity_diff %>%
        diff_vector_2_df(source = "char_height")

      xcal = aa_formula$xcalibur$ratios$nap_intensity_diff %>%
        diff_vector_2_df(source = "xcalibur")

      msnbase = aa_formula$msnbase$ratios$nap_intensity_diff %>%
        diff_vector_2_df(source = "msnbase")

      all_data = rbind(corrected_scan10_max,
                       corrected_char,
                       raw_scan10_max,
                       raw_char,
                       xcal,
                       msnbase)

      all_data = dplyr::left_join(all_data, n_missing, by = "ratio")

      do_differences = list(xcal_raw_char = c("xcalibur", "char_height"),
                            xcal_raw_scan10 = c("xcalibur", "scan_max_raw"),
                            xcal_corrected_char = c("xcalibur",
                                                    "char_corrected"),
                            xcal_corrected_scan10 = c("xcalibur",
                                                      "scan_max_corrected"),
                            char_raw_corrected = c("char_height", "char_corrected"),
                            scan_raw_corrected = c("scan_max_raw", "scan_max_corrected"))

      method_differences = purrr::imap_dfr(do_differences,
                                           function(in_diffs, diff_name){
                                             d_denom = all_data %>%
                                               dplyr::filter(source %in% in_diffs[2])
                                             d_num = all_data %>%
                                               dplyr::filter(source %in% in_diffs[1])

                                             out_diff = data.frame(ratio = d_denom$ratio,
                                                                   diff = d_num$diff - d_denom$diff,
                                                                   source = diff_name)
                                             out_diff
                                           })

      method_differences = dplyr::left_join(method_differences, n_missing, by = "ratio")

      list(data_ratios = all_data,
           method_differences = method_differences)

    })

  })


}

find_large_diffs = function(aa_diffs){
  check_largish = purrr::imap(aa_diffs, function(in_diff, diff_id){
    #message(diff_id)
    purrr::imap(in_diff, function(check_diff, in_id){
      if (is.null(check_diff)) {
        return(NULL)
      }
      #message(in_id)
      check_diff$method_differences %>%
        dplyr::arrange(dplyr::desc(diff)) %>%
        dplyr::slice_head(n = 1)
    })
  })
  check_largish
}

large_figure = function(difference_data, nap_height_data, xcalibur_data, aa_id = "Threonine", formula_id = "C9H17N1Na1O5.15N.0", use_data = "raw"){
  # difference_data = tar_read(aa_compared_filtersd_1ecf)
  # nap_height_data = tar_read(nap_height_1ecf)
  # aa_id = "Threonine"
  # formula_id = "C9H17N1Na1O5.15N.0"
  diff_formula = difference_data[[aa_id]][[formula_id]]
  nap_height_formula = nap_height_data[[aa_id]][[formula_id]]

  scan_level = nap_height_formula$scanlevel$ratios$Log10Height$nap_intensity
  n_na = data.frame(ratio = seq(1, ncol(scan_level)),
                    n_missing = apply(scan_level, 2, function(.x){
                      sum(is.na(.x))
                    })) %>%
    dplyr::arrange(n_missing) %>%
    dplyr::mutate(miss_order = seq(1, nrow(.)))

  raw_ratios_char = nap_height_formula$characterized$ratios$Log10Height$nap_intensity_diff %>%
      diff_vector_2_df(source = "Log10Height")
  raw_ratios_xcal = nap_height_formula$xcalibur$ratios$nap_intensity_diff %>%
      diff_vector_2_df(source = "xcalibur")
  raw_ratios_diff = dplyr::left_join(raw_ratios_char,
                                     raw_ratios_xcal, by = "ratio", suffix = c(".char", ".xcal")) %>%
    dplyr::mutate(xcal_char = diff.xcal - diff.char)

  raw_ratios_diff = dplyr::left_join(raw_ratios_diff, n_na, by = "ratio")

  raw_ratios_long = rbind(raw_ratios_char,
                          raw_ratios_xcal) %>%
    dplyr::left_join(., n_na, by = "ratio")

  threonine_rawplot = ggplot(raw_ratios_long, aes(x = ratio, y = diff, color = source, size = n_missing)) +
    geom_point() +
    theme(legend.position = c(0.7, 0.8)) +
    labs(x = "Pairwise Comparison",
         y = "Peak-Peak NAP - Height Difference")

  threonine_diffplot = ggplot(raw_ratios_diff, aes(x = n_missing, y = xcal_char)) +
    geom_point(size = 2) +
    labs(x = "Number of Missing Scans",
         y = "Xcalibur - Characterized NAP - Height Difference")




}
