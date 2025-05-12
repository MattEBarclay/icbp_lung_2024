
# Wilson score CI function ------------------------------------------------

wilson_lower <- function(n, n_trt) {
  lower_ci <- (1/(1+(qnorm(0.975)^2)/n))*((n_trt/n)+(qnorm(0.975)^2)/(2*n)) - (qnorm(0.975)/(1+((qnorm(0.975)^2)/n)))*sqrt((n_trt/n)*(1-n_trt/n)/n + (qnorm(0.975)^2)/(4*(n^2)))
}
wilson_upper <- function(n, n_trt) {
  lower_ci <- (1/(1+(qnorm(0.975)^2)/n))*((n_trt/n)+(qnorm(0.975)^2)/(2*n)) + (qnorm(0.975)/(1+((qnorm(0.975)^2)/n)))*sqrt((n_trt/n)*(1-n_trt/n)/n + (qnorm(0.975)^2)/(4*(n^2)))
}


# Source the data loading file-------------------------------------------------------------------------
source("../icbp_m9_trt_data/meta0_loading_data.R")

# remove small Ns from data files, result files of paper --------------------

# when N < 10,
# - suppress n
# - suppress N
# Round other ns and Ns in same breakdown to nearest 10

# if total less than 10
trt_unadj_n <- trt_unadj_n |>
  mutate(flag = n < 10 & n != 0) |>
  group_by(trt, jurisdiction, stage, cancer, variable) |>
  mutate(flag_max = max(flag, na.rm = TRUE))

trt_unadj_n <- trt_unadj_n |>
  mutate(n     = ifelse(flag_max, 10*round(n    /10,0), n    )) |>
  mutate(n     = ifelse(flag, NA, n)) |>
  mutate(n_trt = ifelse(flag_max, 10*round(n_trt/10,0), n_trt)) |>
  mutate(n_trt = ifelse(flag, NA, n_trt)) 

trt_unadj_n <- trt_unadj_n |>
  mutate(n_hi     = ifelse(flag_max, n    +4, n)) |>
  mutate(n_lo     = ifelse(flag_max, n    -5, n)) |>
  mutate(n_trt_hi = ifelse(flag_max, n_trt+4, n_trt)) |>
  mutate(n_trt_lo = ifelse(flag_max, n_trt-5, n_trt))

trt_unadj_n |> filter(cancer == "Lung") |> filter(flag_max == 1) 
# no impact for lung

# if treated less than 10
trt_unadj_n <- trt_unadj_n |>
  mutate(flag = n_trt < 10 & n_trt != 0) |>
  group_by(trt, jurisdiction, stage, cancer, variable) |>
  mutate(flag_max = max(flag, na.rm = TRUE))

trt_unadj_n <- trt_unadj_n |>
  mutate(n_trt = ifelse(flag_max, 10*round(n_trt/10,0), n_trt)) |>
  mutate(n_trt = ifelse(flag, NA, n_trt)) 

trt_unadj_n <- trt_unadj_n |>
  mutate(n_trt_hi = ifelse(flag_max, n_trt+4, n_trt)) |>
  mutate(n_trt_lo = ifelse(flag_max, n_trt-5, n_trt))

trt_unadj_n |> filter(cancer == "Lung") |> filter(flag_max == 1) |> select(jurisdiction) |> distinct() |> print(n = 100)
# small impact for some jurisdictions

# recalculate CIs
trt_unadj_n <- trt_unadj_n |>
  mutate(
    prop = n_trt/n,
    lower = wilson_lower(n = n_hi, n_trt = n_trt_lo),
    upper = wilson_upper(n = n_lo, n_trt = n_trt_hi)
  )

# clean up
trt_unadj_n <- trt_unadj_n |>
  select(-flag, -flag_max, -n_hi, -n_lo, -n_trt_hi, -n_trt_lo) |>
  group_by()

# Save the necessary/useful data files produced by the data loadin --------
trt_q_perc |> filter(variable == "tumour_topography_group", variable_value == "Lung") |> saveRDS(file = "Rdata/trt_q_perc.RDS")
trt_unadj_n |> filter(cancer == "Lung") |> saveRDS(file = "Rdata/trt_unadj_n.RDS")

save(
  country_order, 
  country_order_chm, 
  country_order_reml, 
  country_order_reml2,
  site,
  site_chm,
  term_list,
  file = "Rdata/useful_vectors.Rdata"
  )

rm(list = ls())
