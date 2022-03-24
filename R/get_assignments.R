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

find_aa_assignments = function(assign_data, aa_formulas){
  scored_assignments = smirfeTools::score_filter_assignments(assign_data, filter_conditions = e_value <= 0.5)
  within_emfs = smirfeTools:::get_sample_emfs(scored_assignments$assignments, assign_data$sample, evalue_cutoff = 0.5, use_corroborating = TRUE)

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

  split_aa = split(all_aa, all_aa$AA)
  list(aa = split_aa,
       aa_formulas = aa_formulas,
       scored = scored_assignments)
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
  comp_df
}

process_imf = function(imf_string){
  element_pieces = stringr::str_extract_all(imf_string, '(\\d+[A-Z][a-z]?)(\\d+)')[[1]]
  isotopes = unlist(stringr::str_extract_all(element_pieces, '\\d+[A-Z][a-z]?'))
  counts = stringr::str_remove_all(element_pieces, '\\d+[A-Z][a-z]?')
  data.frame(isotopes = isotopes, count = counts)
}

aa_motivation = function(filtersd, xcalibur, msnbase){
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
    # Glutamine is buggy here and not being split up correctly
    message(aa_id)

    imf_n14n15 = purrr::map_df(use_aa$complete_IMF, function(.x){
      #message(.x)
      iso_data = process_imf(.x)
      if (any(grepl("15N", iso_data$isotopes))) {
        n15_count = as.numeric(iso_data[iso_data$isotopes %in% "15N", "count"])
      } else {
        n15_count = 0
      }
      if (any(grepl("14N", iso_data$isotopes))) {
        n14_count = as.numeric(iso_data[iso_data$isotopes %in% "14N", "count"])
      } else {
        n14_count = 0
      }
      out_data = data.frame(N.15 = n15_count,
                            N.14 = n14_count)
      out_data
    })
    imf_n14n15 = imf_n14n15 %>%
      dplyr::mutate(pattern = paste0(N.14, "_", N.15))
    imf_n14n15$complete_IMF = use_aa$complete_IMF
    use_aa = dplyr::left_join(use_aa, imf_n14n15, by = "complete_IMF")
    use_aa = use_aa %>%
      dplyr::mutate(complete_label = paste0(complete_EMF, ".", lbl.type, ".", lbl.count))


    n15n14_rle = rle(use_aa$pattern)
    label_rle = rle(use_aa$complete_label)

    if (!isTRUE(all.equal(n15n14_rle$lengths, label_rle$lengths))) {
      message("N14 / N15 pattern didn't match labeling!")
    }

    aa_data = dplyr::left_join(dplyr::select(use_aa, -PeakID, -Sample, -ObservedMZ), filtersd$scored$data, by = "Sample_Peak")

    if (aa_id %in% c("Glutamine")) {
      aa_complete = split(aa_data, aa_data$pattern)
    } else {
      aa_complete = split(aa_data, aa_data$complete_label)
    }


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
                             Sample_Peak = Sample_Peak,
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

