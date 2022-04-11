#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param coefficients_data
#' @return
#' @author rmflight
#' @export
coefficients_analysis <- function(coefficients_data) {

  n_freq = purrr::map_dbl(coefficients_data, function(.x){
    length(.x$frequency$bad)
  })

  example_freq = coefficients_data[[(which(n_freq == 7)[1])]]

  df_freq = data.frame(coefficient = example_freq$frequency$coefficients,
                       sample = gsub(".zip$", "", basename(example_freq$zip)))
  df_freq = df_freq %>%
    dplyr::mutate(outlier = coefficient %in% example_freq$frequency$bad)

  n_norm = purrr::map_dbl(coefficients_data, function(.x){
    length(.x$normalization)
  })

  list(n_coef_out = n_freq,
       example = df_freq)

}
