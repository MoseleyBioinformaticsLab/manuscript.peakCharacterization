#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param data_97lipid
#' @return
#' @author rmflight
#' @export
compare_noise_cutoff <- function(data_97lipid) {

  noperc = data_97lipid$clone(deep = TRUE)
  noperc$zip_ms$peak_finder$peak_regions$peak_regions =
    reduce_removing_zero(noperc$peak_finder$peak_regions$sliding_regions,
                         noperc$peak_finder$peak_regions$frequency_point_regions$frequency,
                         noperc$peak_finder$quantile_multiplier,
                         noperc$peak_finder$n_point_region)

  perc99 = data_97lipid$clone(deep = TRUE)
  perc99$zip_ms$peak_finder$reduce_sliding_regions()
  return(c(noperc = length(noperc$zip_ms$peak_finder$peak_regions$peak_regions),
           perc99 = length(perc99$zip_ms$peak_finder$peak_regions$peak_regions)))
}
