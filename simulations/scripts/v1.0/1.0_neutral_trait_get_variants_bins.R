library(data.table)

setwd("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data")

dat <- fread("1.0_neutral_out.variant_info.csv")

breaks <- c(
  0, 25, 50, 100, 250, 500,
  1000, 2000, 5000, 20000, Inf
)

dat[, BIN := cut(age, breaks = breaks)]

bin_names <- c(
  "(0,25]"        = "age_0_25",
  "(25,50]"       = "age_25_50",
  "(50,100]"      = "age_50_100",
  "(100,250]"     = "age_100_250",
  "(250,500]"     = "age_250_500",
  "(500,1e+03]"   = "age_500_1000",
  "(1e+03,2e+03]" = "age_1000_2000",
  "(2e+03,5e+03]" = "age_2000_5000",
  "(5e+03,2e+04]" = "age_5000_20000",
  "(2e+04,Inf]"   = "age_20000_inf"
)

dat[, {
  filename <- paste0(
    "1.0_neutral_out.variant_info_",
    bin_names[as.character(.BY$BIN)],
    ".csv"
  )

  fwrite(.SD[, .(site_id)], filename)

  NULL
}, by = BIN]