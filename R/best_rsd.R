find_best_rsd = function(rsd_values){
  other_frame = readr::read_delim(
    'processed   | sample
     singlenorm  | 2ecf
     singlenorm  | 97lipid',
    delim = '|', trim_ws = TRUE)

  order_frame = tidyr::expand_grid(
    sample = c("1ecf", "2ecf", "49lipid", "97lipid"),
    processed = c("filtersd",
                   "doublenorm",
                   "singlenorm_int",
                   "singlenorm",
                   "perc99_nonorm",
                   "noperc_nonorm",
                   "msnbase_only")
  )

  rsd_values2 = dplyr::left_join(order_frame, rsd_values, by = c("sample", "processed"))
  filtersd_rsd = which(rsd_values2$processed %in% "filtersd")
  other_loc = c(11, 25)

  rsd_ft = rsd_values2 %>%
    flextable() %>%
    hline(i = filtersd_rsd[2:(length(filtersd_rsd))] - 1) %>%
    colformat_double(digits = 2) %>%
    bold(i = c(other_loc, filtersd_rsd), j = seq(3, 6), part = "body")
  rsd_ft

}
