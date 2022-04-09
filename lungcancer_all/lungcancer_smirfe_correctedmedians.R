# cd smirfe_related/smirfe_code/SMIRFE_assigner
# python3 -m pipenv shell
# python3 -m AssignMany.py 12 "../../smirfe_dbs/none_nosulfur_1600.db FILE '_assigned.json' '[]' '[\"H\", \"Na\", \"K\", \"NH4\"]' '[1]'" /mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2022-04-02/*.zip
# Database limits: C: 130, N: 7, O: 28, S: 0, P: 3, H: 230, max M/Z 1605, adducts: K, Na, NH4, NAP cutoff: 0.4
library(FTMS.peakCharacterization)
library(knitrProgressBar)

zip_files <- dir("/mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2022-04-02",
                      full.names = TRUE, pattern = ".zip$")
pb <- progress_estimated(length(zip_files))

corrected_medians <- purrr::map_df(zip_files, function(in_file){
  zip_ms = ZipMS$new(in_file)
  zip_ms$load_peak_finder()

  out_df = data.frame(sample = zip_ms$id,
             median = median(10^zip_ms$peak_finder$peak_regions$peak_data$CorrectedLog10Height))
  zip_ms$cleanup()
  update_progress(pb)
  out_df
})

saveRDS(corrected_medians, file = "data/data_output/lung_data/corrected_medians_2022-04-02.rds")
