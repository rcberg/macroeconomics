library(tidyverse)
library(readxl)
# these will help a LOT. all available as of Nov 2025.
# reference text 1: data/raw/input-output-tables/reference_text_1.pdf
# reference text 2: data/raw/input-output-tables/reference_text_2.pdf
# reference text 3: data/raw/input-output-tables/reference_text_3.pdf

# these will help too
exclude_codes = c( "4200ID", "S00401","S00402","S00300","S00900","T005","V00100","V00200","V00300","T006","T008","T005","V00100","T00OTOP","V00300","VABAS","T018","T00TOP","T00SUB","VAPRO","T017" )

commodity_dictionary = 
  read_excel(
    "data/raw/input-output-tables/detail-level/use.xlsx",
    sheet = '2017',
  skip = 5 ) |> 
  select(`Code`,`Commodity Description`)

dictionary = 
  read_excel(
    "data/raw/input-output-tables/detail-level/supply.xlsx",
    sheet = '2017',
    skip = 4,
  n_max = 1) |> 
    select(!`Commodities/Industries`:...2) |> 
    t() %>%
    tibble(
      `Industry Description` = rownames(.),
      Code = .
    ) |> 
  full_join( commodity_dictionary ) |> 
  select(Code, `Industry Description`, `Commodity Description`)

# uses IO data from BEA for 2017, see: https://www.bea.gov/itable/input-output
use_df = 
  read_excel(
    "data/raw/input-output-tables/detail-level/use.xlsx",
    sheet = '2017',
    skip = 5
  ) |> 
  mutate( across(everything(), ~replace_na(.,0))) 

use_matrix_total = 
  filter( use_df, !grepl("V|T",Code) ) |> 
  select( `1111A0`:`S00203`) |>
  mutate( across(everything(), ~replace_na(.,0))) |> 
  as.matrix() |> 
  unname()

use_matrix = 
  filter( use_df, !(Code %in% exclude_codes) ) |> 
  select( `1111A0`:`S00203`, -`4200ID`) |>
  mutate( across(everything(), ~replace_na(.,0))) |> 
  as.matrix() |> 
  unname()

supply_df = 
  read_excel(
    "data/raw/input-output-tables/detail-level/supply.xlsx",
    sheet = '2017',
    skip = 5
  )|> 
  mutate( across(everything(), ~replace_na(.,0))) 

imports_col = 
  supply_df$MCIF[!(supply_df$Code %in% exclude_codes)] + 
  supply_df$MADJ[!(supply_df$Code %in% exclude_codes)]

make_matrix_total = 
  filter(supply_df, !Code %in% c("4200ID","T017") ) |> 
  select( `1111A0`:`S00203`,-`4200ID`) |>
  as.matrix() |> 
  t()

make_matrix = 
  filter(supply_df, !Code %in% exclude_codes ) |> 
  select( `1111A0`:`S00203`,-`4200ID`) |>
  as.matrix() |> 
  t()

industry_output_col = 
  use_df |>
  filter(Code == "T018") |> 
  select(`1111A0`:`S00203`, -`4200ID`) |> 
  as.numeric()

#scrap_col = 
#  filter( supply_df, Code == "S00401" ) |> 
#  select(`1111A0`:`S00203`, -`4200ID`) |> 
#  as.numeric()
#
#industry_output_col = 
#  industry_output_total - scrap_col
#
#nonscrap_ratio = industry_output_col/industry_output_total
#
#nonscrap_ratio_hat = diag(nonscrap_ratio)
#nonscrap_ratio_hat.inv = solve(nonscrap_ratio_hat)

g_hat = diag(industry_output_col)
g_hat.inv = solve(g_hat)

commodity_output_col = 
  supply_df$T007[!(supply_df$Code %in% exclude_codes)]

q_hat = 
  diag(commodity_output_col)
q_hat.inv =
  solve(q_hat)

direct_input_coef_matrix = 
  use_matrix%*%g_hat.inv

# BEA adjusts the supply table for scrap already, no need for us to adjust this make matrix
mkt_share_matrix = 
  make_matrix%*%q_hat.inv

direct_coef_matrix = 
  direct_input_coef_matrix%*%mkt_share_matrix

com_total_req_matrix_uninv = 
  diag(nrow(direct_coef_matrix)) - direct_coef_matrix
com_total_req_matrix = # commodity x commodity
  solve(com_total_req_matrix_uninv)

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
  mkt_share_matrix%*%com_total_req_matrix

ind_total_req_matrix = # industry x industry
  mkt_share_matrix%*%direct_input_coef_matrix

# validating the total requirements tables
# see last paragraph on p.'12-14' of reference text 1
final_uses_purchaser_value = # same as adding up final uses, BUt it is in puchaser prices
  use_df$T019[!(use_df$Code %in% exclude_codes)] - 
  use_df$T001[!(use_df$Code %in% exclude_codes)]

# a hard lesson learned; BEA chooses purchasers prices for us in this "redefined" use table
# headaches caused by not adjusting for this

purch_price_adj = # allows us to convert to the producer prices used in make_matrix
  (supply_df$T014[!(supply_df$Code %in% exclude_codes)] + 
    supply_df$T015[!(supply_df$Code %in% exclude_codes)] )

final_uses_col = # if we don't do this the validation won't be correct
  final_uses_purchaser_value -
  purch_price_adj -
  imports_col

tot_req_validate_df = 
  tibble(
    check_column = 
      c(
      as.numeric(com_total_req_matrix%*%final_uses_col ), 
      as.numeric(ind_com_total_req_matrix%*%final_uses_col )
    ),
    actual_column = c(commodity_output_col, industry_output_col),
    error_rate_vs_actual = (check_column/actual_column) - 1,
    error_z = abs(check_column - actual_column),
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

# just an error plot
#ggplot( tot_req_validate_df, aes( x = error_rate_vs_actual, fill = output_type) ) + 
#  geom_histogram( bins = 20, color = NA, position = 'dodge') + 
#  labs( x = "Error rate (% from actual)", y = "Number of output categories", fill = "Output type") +
#  scale_x_continuous(labels = scales::percent_format()) +
#  theme_minimal()

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
#  select(`1111A0`:`S00203`,-`4200ID`) |> 
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
#  select(`1111A0`:`S00203`,-`4200ID`) |> 
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
#  filter( Code != "T009", Code != "4200ID" ) |> 
#  select(!1:2, -`4200ID`) |> 
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
#text(0.4, 5e4, bquote(sigma == .(round(sd(as.vector(com_error_matrix)), 5 ))))
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
#text(0.95, 60000, bquote(sigma == .(round(sd(as.vector(ind_error_matrix)), 5 ))))
#
#tibble(
#  computed = as.vector(ind_total_req_matrix),
#  bea = as.vector(ind_comp_total_req_matrix)
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
