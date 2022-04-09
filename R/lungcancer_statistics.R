create_sample_info_df = function(sample_file,
                                 bad_samples,
                                 patient_file){
  # tar_load(sample_file)
  # tar_load(bad_samples)
  # tar_load(patient_file)
  sample_info = read.table(sample_file, header = TRUE, sep = "\t")

  patient_subtype = read.table(patient_file, sep = ",", header = TRUE,
                               stringsAsFactors = FALSE)
  patient_subtype$Patient = as.character(patient_subtype$Patient)

  bad_str = paste0(bad_samples, collapse = "|")
  sample_info = sample_info %>%
    dplyr::filter(!grepl(bad_str, sample))
  sample_info = dplyr::mutate(sample_info, patient = stringr::str_extract(sample, "(uk)?[[:digit:]]*"))
  dup_sample = sample_info$sample[duplicated(sample_info$sample)]

  sample_info = sample_info %>%
    dplyr::filter(!(sample %in% dup_sample))

  sample_info = dplyr::left_join(sample_info, patient_subtype, by = c("patient" = "Patient"))

  sample_info[sample_info$disease == "non-cancer", "Cancer.Type"] = "non-cancer"

  granuloma_patients = dplyr::filter(sample_info, (Cancer.Type %in% "Granuloma")) %>%
    dplyr::pull(patient)

  sample_info = dplyr::filter(sample_info, !(patient %in% granuloma_patients))

  sample_info
}

extract_scancentric_imfs = function(emf_file, which = "raw"){
  emf_data = readRDS(emf_file)

  if (which %in% "raw") {
    imf_height = smirfeTools::extract_imf_emf_data(emf_data, by = "IMF", intensity = Height)
  } else {
    imf_height = smirfeTools::extract_imf_emf_data(emf_data, by = "IMF", intensity = CorrectedLog10Height)
    imf_height$intensity = 10^imf_height$intensity
  }

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
  other_imfs$location = other_locations
  other_imfs
}


get_scancentric_medians = function(scancentric_json){
  scancentric_data = readRDS(scancentric_json)
  medians = purrr::imap_dfr(scancentric_data, function(.x, .y){
    data.frame(sample = .y, median = .x$peak$other_info$median_intensity,
               stringsAsFactors = FALSE)
  })
  medians
}

get_other_medians = function(other_file){
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
  medians = purrr::map_df(other_data, function(.x){
    data.frame(sample = .x$sample_id[1],
               median = median(.x$intensity))
  })
  medians
}

scancentric_qcqa = function(scancentric_imfs, scancentric_medians, sample_info){
  # tar_load(scancentric_imfs)
  # tar_load(scancentric_medians)
  # tar_load(sample_info)
  pos_intensity = scancentric_imfs$intensity
  use_samples = data.frame(sample = intersect(sample_info$sample, colnames(pos_intensity)))
  sample_info = left_join(use_samples, sample_info, by = "sample")

  norm_factors = scancentric_medians$median
  names(norm_factors) = scancentric_medians$sample

  normalized_intensity = metabolomicsUtilities::mu_apply_normalization(pos_intensity, norm_factors)

  all_intensities = data.frame(intensity = as.vector(pos_intensity))
  threshold_value = metabolomicsUtilities::mu_calculate_threshold(all_intensities) / max(norm_factors)

  normalized_intensity[is.na(normalized_intensity)] = threshold_value

  normalized_intensity = normalized_intensity[, sample_info$sample]

  normalized_intensity_disease = t(keep_non_zero_percentage(t(normalized_intensity), sample_info$disease, 0.25, zero_value = threshold_value))

  log_intensity = log2(normalized_intensity_disease)
  pca_decomp = prcomp(t(log_intensity), center = TRUE, scale. = FALSE)
  pca_scores = cbind(as.data.frame(pca_decomp$x), sample_info)

  pca_score_contrib = visqc_score_contributions(pca_decomp$x)

  sample_cor_kt = ICIKendallTau::ici_kendalltau(t(log_intensity), global_na = c(NA, Inf, log2(threshold_value)),
                                 perspective = "global", scale_max = TRUE, diag_good = TRUE)
  sample_cor = sample_cor_kt$cor
  sample_cor[is.na(sample_cor)] = 0
  rownames(sample_info) = sample_info$sample
  sample_order = similarity_reorderbyclass(sample_cor, sample_info[, "disease", drop = FALSE], transform = "sub_1")
  sample_medcor = median_correlations(sample_cor, sample_info$disease)
  sample_outlier = visualizationQualityControl::determine_outliers(median_correlations = sample_medcor)

  list(keep_samples = sample_outlier %>%
         dplyr::filter(!outlier) %>%
           dplyr::pull(sample_id),
       outlier_data = sample_outlier,
       correlation = sample_cor_kt,
       pca = list(scores = pca_scores,
                  contribution = pca_score_contrib))
}

run_imf_models = function(imf_data,
                          median_data,
                          scancentric_qcqa,
                          sample_info,
                          id){
  # imf_data = tar_read(scancentric_imfs_raw)
  # median_data = tar_read(scancentric_medians)
  # scancentric_qcqa = tar_read(qcqa)
  # tar_load(sample_info)

  not_outliers = scancentric_qcqa$outlier_data %>%
    dplyr::filter(!outlier) %>%
    dplyr::pull(sample_id)
  pos_intensity = imf_data$intensity[, not_outliers]
  use_samples = data.frame(sample = intersect(sample_info$sample, colnames(pos_intensity)))
  sample_info = left_join(use_samples, sample_info, by = "sample")

  norm_factors = median_data$median
  names(norm_factors) = median_data$sample

  normalized_intensity = metabolomicsUtilities::mu_apply_normalization(pos_intensity, norm_factors)
  missing_intensity = normalized_intensity
  missing_intensity[is.na(missing_intensity)] = 0
  missing_keep = t(keep_non_zero_percentage(t(missing_intensity), sample_info$disease, 0.5, zero_value = 0, all = TRUE))

  disease_status = factor(sample_info$disease)
  instrument_status = factor(sample_info$instrument)

  missing_log = log(missing_keep)
  missing_log[is.infinite(missing_log)] = NA

  t_res = genefilter::rowttests(missing_log, disease_status, na.rm = TRUE)
  n_miss = apply(missing_log, 1, function(.x){sum(is.na(.x))})
  t_res$imf = rownames(t_res)
  t_res$n_missing = n_miss
  t_res$adj.p.value = p.adjust(t_res$p.value, method = "BH")
  t_res$id = id
  t_res

}



compare_ttests = function(compare_list,
                          reference = "raw"){
  ref_df = compare_list[[reference]] %>%
    dplyr::mutate(log_p = -1 * log10(p.value),
                  log_padjust = -1 * log10(adj.p.value))

  compare_results = purrr::map_df(compare_list, function(in_df){
    if (in_df$id[1] %in% reference) {
      return(NULL)
    }
    test_df = in_df %>%
      dplyr::mutate(log_p = -1 * log10(p.value),
                    log_padjust = -1 * log10(adj.p.value))
    diff_df = dplyr::left_join(test_df, ref_df, by = "imf",
                               suffix = c(".test", ".ref"))

    res_df = diff_df %>%
      dplyr::mutate(p_diff = log_p.test - log_p.ref,
                       padjust_diff = log_padjust.test - log_padjust.ref,
                       nmiss_diff = n_missing.test - n_missing.ref)
    res_df
  })

  compare_results
}

binomial_compare = function(lung_compare, p_use = "p_diff"){
  split_compare = split(lung_compare[[p_use]], lung_compare$id.test)

  binom_res = purrr::map_df(split_compare, function(in_compare){
    n_positive = sum(in_compare > 0)
    n_negative = sum(in_compare < 0)
    n_total = length(in_compare)
    test_res = broom::tidy(binom.test(n_positive, n_total, p = 0.5, alternative = "two.sided"))

    test_res = test_res %>%
      dplyr::mutate(n_positive = n_positive,
                    n_negative = n_negative,
                    direction = dplyr::case_when(
                      n_positive > n_negative ~ "positive",
                      n_positive < n_negative ~ "negative"
                    ))
  })
  binom_res$comparison = names(split_compare)
  binom_res$p.adjust = p.adjust(binom_res$p.value, method = "bonferroni")
  binom_res
}

ttest_compare_diffs = function(lung_compare, p_use = "p_diff"){
  split_compare = split(lung_compare[[p_use]], lung_compare$id.test)

  ttest_res = purrr::map_df(split_compare, function(in_compare){
    test_res = broom::tidy(t.test(in_compare, alternative = "two.sided"))

  })
  ttest_res$comparison = names(split_compare)
  ttest_res$p.adjust = p.adjust(ttest_res$p.value, method = "bonferroni")
  ttest_res
}
