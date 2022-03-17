# python3 -m AssignMany.py 12 "python3 ./Main.py /mlab/scratch/cesb_data/smirfe_dbs/none_nosulfur_1600.db FILE '_assigned.json' '[]' '[\"H\", \"NH4\", \"Na\", \"K\"]' '[1]'" /mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2020-01-27/*.zip
# python3 -m AssignMany.py 12 "python3 ./Main.py /mlab/scratch/cesb_data/smirfe_dbs/none_1600.db FILE '_assigned_sulfur.json' '[]' '[\"H\", \"NH4\", \"Na\", \"K\"]' '[1]'" /mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2020-01-27/*.zip
# Example Workflow Across 3 Groups
# Database limits: C: 130, N: 7, O: 28, S: 0, P: 3, H: 230, max M/Z 1605, adducts: K, Na, NH4, NAP cutoff: 0.4
library(smirfeTools)
library(ggplot2)
library(dplyr)
library(furrr)
library(metabolomicsUtilities)
library(knitrProgressBar)
theme_set(cowplot::theme_cowplot())
plan(multicore)
set_internal_map(furrr::future_map)
options(future.globals.maxSize = 800 * 1024 ^ 2)

assigned_files <- dir("/mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2022-03-08",
                      full.names = TRUE, pattern = "assigned.json")
pb <- progress_estimated(length(assigned_files))

assigned_data <- purrr::map(assigned_files, function(in_file){
  tmp_assign = read_smirfe_assignment(in_file)
  filter_data = score_filter_assignments(tmp_assign, filter_conditions = e_value <= 0.5)
  update_progress(pb)
  filter_data
})

names(assigned_data) = purrr::map_chr(assigned_data, ~ .x$sample)

all_emfs = unique(unlist(purrr::map(assigned_data, ~ .x$assignments$isotopologue_EMF)))

saveRDS(assigned_data, "lung_matched_tissue_raw_smirfe_assignments_2020-01-27.rds")

cat(all_emfs, sep = "\n", file = "all_emfs_2020-01-27.txt")

zip_files <- dir("/mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2020-01-27",
                 full.names = TRUE, pattern = ".zip$")
all_coefficients = extract_coefficient_data(zip_files)

names(all_coefficients) = purrr::map_chr(all_coefficients, "sample")


coefficient_df = purrr::map_df(all_coefficients, function(in_sample){
  data.frame(sample = in_sample$sample,
             sqrt = in_sample$coefficients$frequency_coefficients[2],
             stringsAsFactors = FALSE)
})

coefficient_df = dplyr::mutate(coefficient_df, cluster =
                                 dplyr::case_when(
                                   sqrt < 29800000 ~ 1,
                                   dplyr::between(sqrt, 29802000, 29802600) ~ 2,
                                   sqrt > 29802600 ~ 3
                                 ))

assigned_data = assigned_data[coefficient_df$sample]

grouped_data = split(assigned_data, coefficient_df$cluster)

freq_sd = sd_by_group(grouped_data)

sd_mode = calculate_mode(freq_sd$sd) * 2

classified_emfs = import_emf_classifications("all_emfs_classified_2020-01-27.json")
lipid_df = weight_lipid_classifications(classified_emfs, lipid_weight = 2, not_lipid_weight = 2)
assigned_data = purrr::map(assigned_data, function(in_data){
  score_filter_assignments(in_data, filter_conditions = e_value <= 0.5,
                           emf_weight = lipid_df)
})

grouped_data2 = split(assigned_data, coefficient_df$cluster)

grouped_mz = grouped_mz_after_freq(grouped_data2, sd_mode)

#ggplot(dplyr::filter(grouped_mz, Value <= 0.1), aes(x = Index, y = Value)) + geom_point() + facet_wrap(~ set, ncol = 1)

#ggplot(dplyr::filter(grouped_mz, Value <= 0.1), aes(x = Index, y = Value, color = as.factor(set))) + geom_point(alpha = 0.5)

grouped_mz = dplyr::filter(grouped_mz, Value <= 0.1)

mz_cutoff = fit_predict_mz_cutoff(grouped_mz)

all_vote = extract_assigned_data(assigned_data, difference_cutoff = mz_cutoff,
                                 difference_measure = "ObservedMZ", progress = TRUE)
saveRDS(all_vote, file = "lung_voted_all_2020-01-27.rds")

textme::textme("Lung is all done!")
