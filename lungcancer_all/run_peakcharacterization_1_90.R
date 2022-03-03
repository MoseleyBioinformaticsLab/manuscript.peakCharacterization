library(methods)
library(FTMS.peakCharacterization)
pkg = utils::packageDescription("FTMS.peakCharacterization")
# actually run this stuff

peak_pick_samples <- read.table(here::here("lungcancer_all", "file_sample_info.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
all_mzml_files <- dir("/mlab/data/archives/FTMS_raw/CESB/mzml_data", recursive = TRUE, full.names = TRUE, pattern = "mzML$")

all_mzml_files <- data.frame(mzml_file = basename(all_mzml_files), mzml_path = all_mzml_files, stringsAsFactors = FALSE, sample = gsub(".mzML", "", basename(all_mzml_files)))

peak_pick_samples <- dplyr::left_join(peak_pick_samples, all_mzml_files, by = "sample")

curr_date = Sys.Date()
use_date = "2022-03-03"
if (curr_date != use_date) {
  warning("current date and the use date dont match!")
}
zip_save <- file.path("/mlab/scratch/cesb_data/zip_files", paste0("lung_matched_tissue-", use_date))
if (!dir.exists(zip_save)) {
  dir.create(zip_save)
}

existing_files <- dir(zip_save, pattern = ".zip", full.names = TRUE)
existing_samples <- gsub(".zip", "", basename(existing_files))

if (!file.exists(file.path(zip_save, "pkg.rds"))) {
  saveRDS(pkg, file.path(zip_save, "pkg.rds"))
} else {
  old_pkg = readRDS(file.path(zip_save, "pkg.rds"))
  if (pkg != old_pkg) {
    warning("saved pkg and current pkg don't match!")
  }
}

peak_pick_samples <- peak_pick_samples[!(peak_pick_samples$sample %in% existing_samples), ]

n_sample <- floor(nrow(peak_pick_samples) / 2)
use_files <- peak_pick_samples$mzml_path[1:n_sample]

library(furrr)
plan(multiprocess)
set_internal_map(furrr::future_map)

json_files = gsub(".mzML", ".json", use_files)
zip_results = run_mzml_list(use_files, json_files = json_files, save_loc = zip_save)

