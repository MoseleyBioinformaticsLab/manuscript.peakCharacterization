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

calculate_log_ratio_differences = function(nap_df){
  nap_df = nap_df %>%
    dplyr::arrange(dplyr::desc(NAP))

  n_comp = nrow(nap_df) * (nrow(nap_df) - 1) / 2
  comp_df = data.frame(i = rep(NA, n_comp),
                       j = rep(NA, n_comp),
                       nap_ratio = rep(NA, n_comp),
                       intensity_ratio = rep(NA, n_comp))
  comp_count = 1
  for (i in seq_len(nrow(nap_df) - 1)) {

    for (j in seq(i + 1, nrow(nap_df))) {
      #message(paste0(i, "_", j))
      comp_df[comp_count, "nap_ratio"] = log(nap_df$NAP[i]) - log(nap_df$NAP[j])
      comp_df[comp_count, "intensity_ratio"] = log(nap_df$intensity[i]) - log(nap_df$intensity[j])
      comp_df[comp_count, "i"] = nap_df$Sample_Peak[i]
      comp_df[comp_count, "j"] = nap_df$Sample_Peak[j]
      comp_count = comp_count + 1
    }

  }
  comp_df = comp_df %>%
    dplyr::mutate(nap_intensity_diff = nap_ratio - intensity_ratio)
  p1 = get_number_sample_peak(comp_df$i)
  p2 = get_number_sample_peak(comp_df$j)

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

aa_height_nap_scanlevel = function(aa_data){
  # aa_data = tar_read(aa_filtersd_1ecf)
  scan_level = aa_data$scored$scan_level
  # make sure to purrr::map here
  aa_peak_data = aa_data$aa$Alanine

  aa_peak_data = aa_peak_data %>%
    dplyr::mutate(complete_label = paste0(complete_EMF, "_", lbl.type, "_", lbl.count))

  split_aa = split(aa_peak_data, aa_peak_data$complete_label)
  split_big = purrr::map_lgl(split_aa, ~ nrow(.x) > 1)
  split_aa = split_aa[split_big]

  # and purrr::map again here
  split_use = split_aa[[1]]
  intensity_use = scan_level$Log10Height[split_use$Sample_Peak]
  intensity_matrix = log(10^as.matrix(dplyr::bind_rows(intensity_use)))

  nap_matrix = matrix(log(split_use$NAP),
                      nrow = nrow(intensity_matrix),
                      ncol = ncol(intensity_matrix),
                      byrow = TRUE)
  colnames(nap_matrix) = colnames(intensity_matrix)

  nap_comparisons = calculate_comparisons_matrix(nap_matrix)
  intensity_comparisons = calculate_comparisons_matrix(intensity_matrix)

  nap_intensity = abs(nap_comparisons$matrix - intensity_comparisons$matrix)
  rownames(nap_intensity) = seq(1, nrow(nap_intensity))

  no_na = apply(nap_intensity, 1, function(.x){
    sum(is.na(.x)) == 0
  })
  nap_intensity = nap_intensity[no_na, ]
  nap_intensity_max = rowMax(abs(nap_intensity))
  nap_intensity_median = rowMedians(abs(nap_intensity))

  order_max = order(nap_intensity_max)
  nap_intensity_bymax = nap_intensity[order_max, ]
  order_median = order(nap_intensity_median)
  nap_intensity = nap_intensity[order_median, ]
}

aa_height_nap_all = function(filtersd, xcalibur, msnbase){
  # filtersd = tar_read(aa_filtersd_1ecf)
  # xcalibur = tar_read(xcalibur_1ecf)
  # msnbase = tar_read(msnbase_1ecf)
  xcalibur_peaks = xcalibur %>%
    dplyr::transmute(mz = ...1,
                  intensity = ...2)
  xcalibur_peaks$Sample_Peak = paste0("xcalibur_", seq_len(nrow(xcalibur_peaks)))

  msnbase_peaks = msnbase$comb
  msnbase_peaks$Sample_Peak = paste0("msnbase_", seq_len(nrow(msnbase_peaks)))

  aa_info = purrr::imap(filtersd$aa, function(use_aa, aa_id){
    # use_aa = filtersd$aa[["Glutamine"]]
    # aa_id = "Glutamine"
    message(aa_id)

    use_aa = use_aa %>%
      dplyr::mutate(complete_label = paste0(complete_EMF, ".", lbl.type, ".", lbl.count))

    aa_complete = split(use_aa, use_aa$complete_label)
    aa_hasenough = purrr::map_lgl(aa_complete, ~ (nrow(.x) > 1))
    aa_complete = aa_complete[aa_hasenough]
    aa_ratios = purrr::imap(aa_complete, function(in_complete, split_id){
      message(split_id)
      use_intensities = c("CorrectedLog10Height",
                      "Height")

      char_level = purrr::map(use_intensities, function(intensity){
        tmp_df = in_complete[, c("Sample_Peak", "NAP", intensity)]
        if (grepl("Log10", intensity)) {
          tmp_df$intensity = 10^tmp_df[[intensity]]
        } else {
          tmp_df$intensity = tmp_df[[intensity]]
        }
        peak_ratios = calculate_log_ratio_differences(tmp_df)
        peak_ratios$intensity_measure = intensity
        peak_ratios$source = "characterized"
        peak_ratios
      })

      # we should be able to do this at the individual scan level
      # as well, but that's more work.
      scan_level = NULL

      complete_theory = in_complete %>%
        dplyr::transmute(NAP = NAP,
                         theoreticalmz = ObservedMZ - mass_error,
                         Sample_Peak = Sample_Peak,
                         isotopologue_IMF = isotopologue_IMF,
                         lowermz = theoreticalmz - (theoreticalmz * 2e-6),
                         highmz = theoreticalmz + (theoreticalmz * 2e-6))

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

      xcal_matched = match_peaks(complete_theory, xcalibur_peaks)

      if (nrow(xcal_matched) > 1) {
        xcal_ratios = calculate_log_ratio_differences(xcal_matched[, c("Sample_Peak", "NAP", "intensity")])
        xcal_ratios$source = "xcalibur"
      } else {
        xcal_ratios = NULL
      }

      msnbase_matched = match_peaks(complete_theory, msnbase_peaks)

      if (nrow(msnbase_matched) > 1) {
        msnbase_ratios = calculate_log_ratio_differences(msnbase_matched[, c("Sample_Peak", "NAP", "intensity")])
        msnbase_ratios$source = "msnbase"
      } else {
        msnbase_ratios = NULL
      }

      list(characterized = list(ratios = char_level,
                                peaks = in_complete),
           xcalibur = list(ratios = xcal_ratios,
                           peaks = xcal_matched),
           msnbase = list(ratios = msnbase_ratios,
                          peaks = msnbase_matched))

    })
  })

}

