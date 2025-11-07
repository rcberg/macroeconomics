library(tidyverse)
library(readxl)
# these will help a LOT. all available as of Nov 2025.
# reference text 1: data/raw/input-output-tables/reference_text_1.pdf
# reference text 2: data/raw/input-output-tables/reference_text_2.pdf
# reference text 3: data/raw/input-output-tables/reference_text_3.pdf

# these will help too
exclude_codes = 
  c( "S00401","S00402","S00300","S00900","T005","V00100","V00200","V00300","T006","T008" )

commodity_dictionary = 
  read_excel(
    "data/raw/input-output-tables/detail-level/use_before_redef_prod_price.xlsx",
    sheet = '2017',
  skip = 5 ) |> 
  select(`Code`,`Commodity Description`)

dictionary = 
  read_excel(
    "data/raw/input-output-tables/detail-level/make.xlsx",
    sheet = '2017',
    skip = 5 ) |>
  dplyr::select(Code, `Industry Description`) |> 
  full_join( commodity_dictionary )

# uses IO data from BEA for 2017, see: https://www.bea.gov/itable/input-output
use_df = 
  read_excel(
    "data/raw/input-output-tables/detail-level/use_before_redef_prod_price.xlsx",
    sheet = '2017',
    skip = 5
  ) |> 
  mutate( across(everything(), ~replace_na(.,0))) 

total_final_uses = 
  use_df$...426[!grepl("V|T",use_df$Code)]

final_uses_col = 
  use_df$...426[!(use_df$Code %in% exclude_codes)]

use_matrix_total = 
  filter( use_df, !grepl("V|T",Code) ) |> 
  select( `1111A0`:`S00203`) |>
  mutate( across(everything(), ~replace_na(.,0))) |> 
  as.matrix() |> 
  unname()

use_matrix = 
  filter( use_df, !(Code %in% exclude_codes) ) |> 
  select( `1111A0`:`S00203`) |>
  mutate( across(everything(), ~replace_na(.,0))) |> 
  as.matrix() |> 
  unname()

make_df = 
  read_excel(
    "data/raw/input-output-tables/detail-level/make.xlsx",
    sheet = '2017',
    skip = 5
  )|> 
  mutate( across(everything(), ~replace_na(.,0))) 

make_matrix_total = 
  filter(make_df, Code != "T007") |> 
  select( `1111A0`:`S00900`) |>
  as.matrix()

make_matrix = 
  make_matrix_total[,-which(colnames(make_matrix_total) %in% c("S00401","S00402","S00300","S00900"))] 

industry_output_total = 
  make_df$T008[make_df$Code != "T007"] 

scrap_col = 
  make_df$S00401[make_df$Code != "T007"]

industry_output_col = 
  industry_output_total - scrap_col

nonscrap_ratio = industry_output_col/industry_output_total

nonscrap_ratio_hat = diag(nonscrap_ratio)
nonscrap_ratio_hat.inv = solve(nonscrap_ratio_hat)

g_hat = diag(industry_output_total)
g_hat.inv = solve(g_hat)

commodity_output_col = 
  use_df$...427[!(use_df$Code %in% exclude_codes)]

q_hat = 
  diag(commodity_output_col)
q_hat.inv =
  solve(q_hat)

direct_input_coef_matrix = 
  use_matrix%*%g_hat.inv

mkt_share_matrix = 
  make_matrix%*%q_hat.inv

transformation_matrix = 
  nonscrap_ratio_hat.inv%*%mkt_share_matrix

direct_coef_matrix = 
  direct_input_coef_matrix%*%transformation_matrix

com_total_req_matrix_uninv = 
  diag(nrow(direct_coef_matrix)) - direct_coef_matrix
com_total_req_matrix = # commodity x commodity
  solve(com_total_req_matrix_uninv) # this seems to work!

#comp_total_req_matrix = 
#  read_excel(
#    "data/raw/input-output-tables/detail-level/com-com_total_requirement.xlsx", 
#    sheet = '2017',skip = 4
#  ) |> 
#  filter( Code != "T010" ) |> 
#  select(`1111A0`:`S00900`) |> 
#  unname() |> 
#  as.matrix()
#
#error_matrix = 
#  comp_total_req_matrix[-c(400:401),-c(400:401)] - total_req_matrix

ind_com_total_req_matrix = # industry x commodity
  transformation_matrix%*%com_total_req_matrix

ind_total_req_matrix = # industry x industry
  transformation_matrix%*%direct_input_coef_matrix

# validating the total requirements tables
# see last paragraph on p.'12-14' of reference text 1

tot_req_validate_df = 
  tibble(
    check_column = 
      c(
      as.numeric(com_total_req_matrix%*%final_uses_col), 
      as.numeric(ind_com_total_req_matrix%*%final_uses_col)
    ),
    actual_column = c(commodity_output_col, industry_output_col),
    error_rate_vs_actual = (check_column/actual_column) - 1,
    output_type = 
      c(
      rep("commodity", length(commodity_output_col)), 
      rep("industry", length(industry_output_col))
     )
  )

ggplot( tot_req_validate_df, aes( x = check_column, y = actual_column) ) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs( x = "Calculated from matrix ($millions)", y = "Actual data ($millions)") +
  scale_x_continuous(labels = scales::dollar_format()) +
  scale_y_continuous(labels = scales::dollar_format()) +
  theme_minimal()

ggplot( tot_req_validate_df, aes( x = error_rate_vs_actual, fill = output_type) ) + 
  geom_histogram( bins = 20, color = NA, position = 'dodge') + 
  labs( x = "Error rate (% from actual)", y = "Number of output categories", fill = "Output type") +
  scale_x_continuous(labels = scales::percent_format()) +
  theme_minimal()

lm( actual_column ~ check_column, tot_req_validate_df ) |> summary()

# plots of errors compared to official BEA total requirements matrix
# feel free to comment these out if not needed

#com_comp_total_req_matrix = 
#  readxl::read_excel(
#    "data/raw/input-output-tables/detail-level/com-com_total_requirement.xlsx",
#    sheet = '2017',
#    skip = 4
#  ) |> 
#  filter( !(Code %in% c(exclude_codes,"T010") )) |> 
#  select(`1111A0`:`S00203`) |> 
#  as.matrix()
#
#com_error_matrix =  # commodity x commodity
#  com_comp_total_req_matrix - com_total_req_matrix
#
#ind_com_comp_total_req_matrix = 
#  readxl::read_excel(
#    "data/raw/input-output-tables/detail-level/ind-com_total_requirement.xlsx",
#    sheet = '2017',
#    skip = 4
#  ) |> 
#  filter( !(Code %in% c(exclude_codes,"T010") )) |> 
#  select(`1111A0`:`S00203`) |> 
#  as.matrix()
#
#ind_com_error_matrix = 
#  ind_com_comp_total_req_matrix - ind_com_total_req_matrix
#
#ind_comp_total_req_matrix = 
#    readxl::read_excel(
#    "data/raw/input-output-tables/detail-level/ind-ind_total_requirement.xlsx",
#    sheet = '2017',
#    skip = 4
#  ) |> 
#  filter( Code != "T009" ) |> 
#  select(!1:2) |> 
#  as.matrix()
#
#ind_error_matrix = 
#  ind_comp_total_req_matrix - ind_total_req_matrix
#
#hist(
#  com_error_matrix,
#  breaks = 30,
#  xlab = "Deviation from BEA",
#  ylab = "Number of matrix entries",
#  main = "Commodity-by-commodity requirement table error"
#)
#text(0.2, 6e4, bquote(sigma == .(round(sd(as.vector(com_error_matrix)), 5 ))))
#
#hist(
#  ind_com_error_matrix,
#  breaks = 30,
#  xlab = "Deviation from BEA",
#  ylab = "Number of matrix entries",
#  main = "Industry-by-commodity requirement table error"
#)
#text(0.35, 50000, bquote(sigma == .(round(sd(as.vector(ind_com_error_matrix)), 5 ))))
#
#hist(
#  ind_error_matrix,
#  breaks = 30,
#  xlab = "Deviation from BEA",
#  ylab = "Number of matrix entries",
#  main = "Industry-by-industry requirement table error"
#)
#text(0.7, 80000, bquote(sigma == .(round(sd(as.vector(ind_error_matrix)), 5 ))))
#
#tibble(
#  computed = as.vector(ind_com_total_req_matrix),
#  bea = as.vector(ind_com_comp_total_req_matrix)
#) |> 
#  ggplot( aes( x = computed, y = bea)) + 
#  geom_point() + 
#  geom_abline(slope = 1, intercept = , linetype = 2) + 
#  theme_minimal()
#
#tibble(
#  computed = as.vector(ind_com_total_req_matrix),
#  bea = as.vector(ind_com_comp_total_req_matrix)
#) |> 
#  lm( formula = bea ~ computed ) |> summary()
