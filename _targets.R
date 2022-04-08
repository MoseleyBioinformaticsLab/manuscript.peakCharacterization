## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

## debugging stuff
#tar_option_set(debug = "scan_height_correlation")
## targets::tar_make(callr_function = NULL)

use_files = dir("data/data_input", pattern = "mzML", full.names = TRUE)
excel_files = dir("data/data_input", pattern = "[ECF|pos].xlsx", full.names = TRUE)

pkg_sha = as.character(utils::packageDescription("FTMS.peakCharacterization")$GithubSHA1)

pkg_tar = tar_target(pkg, pkg_sha)

data_mzml = expand_grid(
  data_function = rlang::syms("reading_scans_tile_windows"),
  mzml = use_files
) %>%
  dplyr::mutate(name = dplyr::case_when(
    grepl("97Cpos", mzml) ~ "97lipid",
    grepl("49Cpos", mzml) ~ "49lipid",
    grepl("1_ECF", mzml) ~ "1ecf",
    grepl("2_ECF", mzml) ~ "2ecf"
  ))

method_data = expand_grid(
  method_function = (c("noperc_nonorm",
                                  "perc99_nonorm",
                                  "singlenorm",
                                  "intsinglenorm",
                                  "doublenorm",
                                  "filtersd")),
  data_names = paste0("data_", data_mzml$name)) %>%
  dplyr::mutate(name = gsub("data_", "", data_names),
                data_sym = rlang::syms(data_names),
                method_sym = rlang::syms(method_function))

data_tar = tar_map(
  values = data_mzml,
  names = "name",
  tar_target(data, data_function(mzml, pkg))
)
method_tar = tar_map(
  values = method_data,
  unlist = FALSE,
  names = c("method_function", "name"),
  tar_target(method, method_sym(data_sym)),
  tar_target(rsd, single_rsd(method)),
  tar_target(zip, write_peaks_for_assignment(method)),
  tar_target(assign, assign_files(zip)),
  tar_target(n_peak, n_final_peak(method)),
  tar_target(emf, get_assignments(assign))
)


normalization_tar = tar_combine(
  normalization_combine,
  method_tar[[1]],
  command = normalization_factors(!!!.x),
  iteration = "list"
)

n_peak_tar = tar_combine(
  n_peak_combine,
  method_tar[[5]],
  command = dplyr::bind_rows(!!!.x)
)

noise_tar = tibble(
  data_name = (c(
    "data_97lipid",
    "data_49lipid",
    "data_1ecf",
    "data_2ecf"
  ))
) %>%
  dplyr::mutate(
    data_syms = rlang::syms(data_name),
    data_id = rename_samples(data_name)
  )

compare_noise_tar = tar_map(
  values = noise_tar,
  names = "data_id",
  tar_target(noise, compare_noise_cutoffs(data_syms))
)

combine_noise_tar = tar_combine(
  noise_combine,
  compare_noise_tar,
  command = dplyr::bind_rows(!!!.x)
)

# rsd_tar = tar_combine(
#   rsd_data,
#   values = ends_with("97lipid")
# )

figures_tar = tar_plan(
  tar_target(frequency_conversion, plot_frequency_conversion(data_97lipid)),
  tar_target(peak_ordering, plot_peak_ordering(method_filtersd_97lipid)),
  tar_target(sliding_regions, plot_sliding_window_density(data_97lipid)),
  tar_target(peak_fit_plot, plot_peak_fitting(method_filtersd_97lipid)),
  tar_target(rsd_plot, plot_rsd_differences(rsd_combine)),

  tar_target(find_sub_region, split_regions(method_perc99_nonorm_97lipid)
  ),
  tar_target(intensity_scan_cor, correlate_scan_height(method_perc99_nonorm_97lipid)),
  tar_target(intensity_scan_cor_plot, correlate_scan_height_graph(intensity_scan_cor)),
  tar_target(split_region_plot,
             plot_region_splitting(find_sub_region)
  ),
  tar_target(compare_normalization,
             normalization_graph(normalization_combine)),
  tar_target(noise_plot,
             create_noise_plot(noise_combine)),
  tar_target(mn_ratios,
             calculate_m_n_ratio(noise_combine)),
  tar_target(motivation_plot,
             motivating_plot(nap_height_1ecf)),
  tar_target(base_manuscript,
             "doc/peakcharacterization_manuscript.Rmd",
             format = "file"),
  tar_render(manuscript_both,
             "doc/peakcharacterization_manuscript.Rmd"),
  tar_target(manuscript_nostyle,
             strip_mdpi_render(base_manuscript)),
  tar_target(manuscript_mdpi,
             strip_headers_render(base_manuscript))
)

tables_tar = tar_plan(
  tar_target(rsd_values, summarize_rsd(rsd_combine))
)

# other_tar = tar_plan(
#   tar_target(
#     nonoise_vs_noise,
#     compare_noise_cutoff(data_97lipid)
#   )
# )

msnbase_mzml = expand_grid(
  data_function = rlang::syms("msnbase_centroid"),
  mzml = use_files
) %>%
  dplyr::mutate(name = dplyr::case_when(
    grepl("97Cpos", mzml) ~ "97lipid",
    grepl("49Cpos", mzml) ~ "49lipid",
    grepl("1_ECF", mzml) ~ "1ecf",
    grepl("2_ECF", mzml) ~ "2ecf"
  ))

msnbase_tar = tar_map(
  values = msnbase_mzml,
  names = "name",
  tar_target(msnbase, data_function(mzml)),
  tar_target(msnbase_scanpeaks, msnbase_match_scan_combined(msnbase)),
  tar_target(rsd_msnbase_scanpeaks, calc_rsd_msnbase(msnbase_scanpeaks))
)

rsd_tar = tar_combine(
  rsd_combine,
  c(method_tar[[2]], msnbase_tar[[3]]),
  command = merge_rsd(!!!.x),
  iteration = "list"
)

hpd_df = tibble(
  method1 = rlang::syms(c("method_filtersd_1ecf",
              "method_filtersd_2ecf",
              "method_filtersd_49lipid",
              "method_filtersd_97lipid")),
  method2 = rlang::syms(c("method_noperc_nonorm_1ecf",
                          "method_noperc_nonorm_2ecf",
                          "method_noperc_nonorm_49lipid",
                          "method_noperc_nonorm_97lipid")),
  method3 = rlang::syms(c("msnbase_1ecf",
                          "msnbase_2ecf",
                          "msnbase_49lipid",
                          "msnbase_97lipid"))
) %>%
  dplyr::mutate(names = rename_samples(method1))

hpd_tar = tar_map(
  values = hpd_df,
  names = "names",
  unlist = FALSE,
  tar_target(hpd, hpds_from_excel(method1, method2, method3, excel_files)),
  #tar_target(hpd_chisq, chisq_hpds(hpd)),
  tar_target(hpd_plots, plot_hpds(hpd))
)

# hpd_width_tar = tar_combine(
#   hpd_width_all,
#   hpd_tar[[4]],
#   iteration = "list",
#   command = combine_hpd_width(!!!.x)
# )
#
# hpd_chisq_tar = tar_combine(
#   hpd_chisq_all,
#   hpd_tar[[2]],
#   command = dplyr::bind_rows(!!!.x)
# )
#
aa_methods = expand_grid(
  assign = "emf",
  method = c("noperc_nonorm",
              "perc99_nonorm",
              "singlenorm",
              "intsinglenorm",
              "doublenorm",
              "filtersd"),
  sample = c("1ecf",
             "2ecf")
) %>%
  dplyr::mutate(in_assign = rlang::syms(paste0(assign, "_", method, "_", sample)),
                name = paste0(method, "_", sample))

formula_tar = tar_target(aa_formula,
                         get_expected_formulas())

aa_tar = tar_map(
  values = aa_methods,
  names = "name",
  tar_target(aa_high, find_assignments(in_assign, aa_formula)),
  tar_target(aa, find_assignments(in_assign, aa_formula, 0.1))
)

lipid_methods = expand_grid(
  assign = "emf",
  method = c("noperc_nonorm",
             "perc99_nonorm",
             "singlenorm",
             "intsinglenorm",
             "doublenorm",
             "filtersd"),
  sample = c("97lipid",
             "49lipid")
) %>%
  dplyr::mutate(in_assign = rlang::syms(paste0(assign, "_", method, "_", sample)),
                name = paste0(method, "_", sample))

lipid_tar = tar_map(
  values = lipid_methods,
  names = "name",
  tar_target(lipid_high, find_assignments(in_assign)),
  tar_target(lipid, find_assignments(in_assign, e_cutoff = 0.1))
)

xcalibur_df = tibble(
  files = excel_files,
  name = rename_samples(excel_files)
)

xcalibur_tar = tar_map(
  values = xcalibur_df,
  names = "name",
  tar_target(xcalibur, get_xcalibur_peaks(files))
)

height_nap_tar = list(
  tar_target(nap_height_1ecf, aa_height_nap_all(aa_filtersd_1ecf, xcalibur_1ecf, msnbase_1ecf)),
  tar_target(nap_height_2ecf, aa_height_nap_all(aa_filtersd_2ecf, xcalibur_2ecf, msnbase_2ecf)),
  tar_target(aa_compared_filtersd_1ecf, compare_peak_ratios(nap_height_1ecf)),
  tar_target(aa_compared_filtersd_2ecf, compare_peak_ratios(nap_height_2ecf)),
  tar_target(aa_large_filtersd_1ecf, find_large_diffs(aa_compared_filtersd_1ecf)),
  tar_target(aa_large_filtersd_2ecf, find_large_diffs(aa_compared_filtersd_2ecf)),
  tar_target(aa_diffs_filtersd_1ecf, ratio_diffs_extracted(aa_compared_filtersd_1ecf)),
  tar_target(aa_diffs_filtersd_2ecf, ratio_diffs_extracted(aa_compared_filtersd_2ecf))
)

unassigned_match_tar = list(
  tar_target(um_peaks_1ecf,
             match_unassigned_peaks(method_filtersd_1ecf, xcalibur_1ecf, msnbase_1ecf)),
  tar_target(um_peaks_2ecf,
             match_unassigned_peaks(method_filtersd_2ecf, xcalibur_2ecf, msnbase_2ecf)),
  tar_target(um_peaks_97lipid,
             match_unassigned_peaks(method_filtersd_97lipid, xcalibur_97lipid, msnbase_97lipid)),
  tar_target(um_peaks_49lipid,
             match_unassigned_peaks(method_filtersd_49lipid, xcalibur_49lipid, msnbase_49lipid))
)

assigned_match_tar = list(
  tar_target(ass_peaks_1ecf,
             match_assigned_peaks(aa_filtersd_1ecf, xcalibur_1ecf, msnbase_1ecf)),
  tar_target(ass_peaks_2ecf,
             match_assigned_peaks(aa_filtersd_2ecf, xcalibur_2ecf, msnbase_2ecf)),
  tar_target(ass_peaks_97lipid,
             match_assigned_peaks(lipid_filtersd_97lipid, xcalibur_97lipid, msnbase_97lipid)),
  tar_target(ass_peaks_49lipid,
             match_assigned_peaks(lipid_filtersd_49lipid, xcalibur_49lipid, msnbase_49lipid))
)

lungcancer_tar = list(
  tar_target(sample_file,
             "data/data_output/lung_data/file_sample_info.txt",
            format = "file"),
  tar_target(patient_file,
             "data/data_output/lung_data/lungcancer_tissue_type_patient_mappings.csv"),
  tar_target(sample_info,
             create_sample_info_df(
               sample_file,
               bad_samples,
               patient_file
             )),
  tar_target(scancentric_json,
             "data/data_output/lung_json_data_2022-04-02.rds",
             format = "file"),
  tar_target(scancentric_medians,
             get_scancentric_medians(scancentric_json)),
  tar_target(emf_file,
             "data/data_output/lung_data/lung_voted_all_2022-04-02.rds",
             format = "file"),
  tar_target(scancentric_imfs_raw,
             extract_scancentric_imfs(
               emf_file
             )),
  tar_target(scancentric_imfs_corrected,
             extract_scancentric_imfs(
               emf_file,
               "corrected"
             )),
  tar_target(msnbase_data,
             "data/data_output/lung_data/lung_msnbase_peaks.rds",
             format = "file"),
  tar_target(msnbase_medians,
             get_other_medians(msnbase_data)),
  tar_target(msnbase_imfs,
             match_imfs(
               msnbase_data,
               scancentric_imfs_raw
               )),
  tar_target(xcalibur_data,
             "data/data_output/lung_data/lung_xcalibur_peaks.rds",
             format = "file"),
  tar_target(xcalibur_medians,
             get_other_medians(xcalibur_data)),
  tar_target(xcalibur_imfs,
             match_imfs(
               xcalibur_data,
               scancentric_imfs_raw
             )),
  tar_target(bad_samples,
             c("151", "154", "158", "70", "77")),
  tar_target(qcqa,
             scancentric_qcqa(
               scancentric_imfs_raw,
               scancentric_medians,
               sample_info)),
  tar_target(lung_raw,
             run_imf_models(
               scancentric_imfs_raw,
               scancentric_medians,
               qcqa,
               sample_info,
               id = "raw"
             )),
  tar_target(lung_corrected,
             run_imf_models(
               scancentric_imfs_corrected,
               scancentric_medians,
               qcqa,
               sample_info,
               id = "corrected"
             )),
  tar_target(lung_msnbase,
             run_imf_models(
               msnbase_imfs,
               msnbase_medians,
               qcqa,
               sample_info,
               id = "msnbase"
             )),
  tar_target(lung_xcalibur,
             run_imf_models(
               xcalibur_imfs,
               xcalibur_medians,
               qcqa,
               sample_info,
               id = "xcalibur"
             )),
  tar_target(lung_compare,
             compare_ttests(
               list(corrected = lung_corrected,
                    raw = lung_raw,
                    msnbase = lung_msnbase,
                    xcalibur = lung_xcalibur),
               reference = "raw")
             ),
  tar_target(lung_binomial,
             binomial_compare(lung_compare))
)

list(pkg_tar,
     data_tar,
     method_tar,
     n_peak_tar,
     compare_noise_tar,
     combine_noise_tar,
     rsd_tar,
     normalization_tar,
     figures_tar,
     tables_tar,
     msnbase_tar,
     hpd_tar,
     formula_tar,
     xcalibur_tar,
     height_nap_tar,
     unassigned_match_tar,
     assigned_match_tar,
     aa_tar,
     lipid_tar,
     lungcancer_tar)
