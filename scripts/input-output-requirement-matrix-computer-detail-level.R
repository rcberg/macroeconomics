library(tidyverse)

# these will help a LOT. all available as of Nov 2025.
# reference text 1: data/raw/input-output-tables/reference_text_1.pdf
# reference text 2: data/raw/input-output-tables/reference_text_2.pdf
# reference text 3: data/raw/input-output-tables/reference_text_3.pdf

# uses IO data from BEA for 2017, see: https://www.bea.gov/itable/input-output
use_df = 
  readxl::read_excel(
    "data/raw/input-output-tables/detail-level/use_before_redef_prod_price.xlsx",
    sheet = '2017',
    skip = 5
  ) |> 
  mutate( across(everything(), ~replace_na(.,0))) 

use_column_dictionary = 
  readxl::read_excel(
    "data/raw/input-output-tables/detail-level/use_before_redef_prod_price.xlsx",
    sheet = '2017' )[4:5,3:404] |> 
  t() |> 
  unname()

unadj_use_matrix = 
  select( use_df, `1111A0`:`S00203`)[1:402,] |> 
  mutate( across(everything(), ~replace_na(.,0))) |> 
  as.matrix() |> 
  unname()

use_matrix = 
  unadj_use_matrix[-c(399,400),]

industry_output_col = 
  select( use_df, `1111A0`:`S00203`)[408,] |> 
  unlist() |> 
  unname()

g_hat = diag(industry_output_col)
g_hat.inv = solve(g_hat)

make_df = 
  readxl::read_excel(
    "data/raw/input-output-tables/detail-level/make.xlsx",
    sheet = '2017',
    skip = 5
  )|> 
  mutate( across(everything(), ~replace_na(.,0))) 

make_column_dictionary = 
  readxl::read_excel(
    "data/raw/input-output-tables/detail-level/make.xlsx",
    sheet = '2017',
    skip = 3,
    n_max = 2 ) |>
  dplyr::select(!1:2) |> 
  t() %>% 
  tibble(
    industry_name = row.names(.),
    industry_code = .
  )

unadj_make_matrix = 
  make_df[1:402,3:404] |>
  unname() |> 
  as.matrix()

make_matrix = 
  unadj_make_matrix[,-c(400:401)] 

commodity_output_col = 
  use_df$...427[c(1:399, 402)] 

q_hat = 
  diag(commodity_output_col)
q_hat.inv =
  solve(q_hat)

direct_req_matrix = 
  use_matrix%*%g_hat.inv

mkt_share_matrix = 
  make_matrix%*%q_hat.inv

direct_coef_matrix = 
  direct_req_matrix%*%mkt_share_matrix

total_req_matrix_uninv = 
  diag(nrow(direct_coef_matrix)) - direct_coef_matrix
total_req_matrix = 
  solve(total_req_matrix_uninv) # this seems to work!

#comp_total_req_matrix = 
#  readxl::read_excel(
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

ind_ind_total_req_matrix = 
  direct_coef_matrix%*%total_req_matrix
