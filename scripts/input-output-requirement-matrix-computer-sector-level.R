library(tidyverse)
# CAUTIONARY NOTE: this 'sector-level' I-O is amazing to practice with, but for many
#                  reasons is BAD to actually use! see reference texts.

# these will help a LOT. all available as of Nov 2025.
# reference text 1: data/raw/input-output-tables/reference_text_1.pdf
# reference text 2: data/raw/input-output-tables/reference_text_2.pdf
# reference text 3: data/raw/input-output-tables/reference_text_3.pdf

# uses IO data from BEA for 2017, see: https://www.bea.gov/itable/input-output
use_df = 
  readxl::read_excel(
    "data/raw/input-output-tables/sector-level/use_before_redef_prod_price.xlsx",
    skip = 8
  ) |> 
  mutate( across(`Agriculture, forestry, fishing, and hunting`:`Total Commodity Output`, ~replace_na(.,0))) 

use_matrix = 
  use_df |> 
  filter( !is.na(Code), !startsWith(Code,"V")) |> 
  select(`Agriculture, forestry, fishing, and hunting`:`Government`) |> 
  as.matrix()

make_df = 
  readxl::read_excel(
    "data/raw/input-output-tables/sector-level/make.xlsx",
    skip = 9
  ) |> 
  mutate( across(`Agriculture, forestry, fishing, and hunting`:`Total Industry Output`, ~replace_na(.,0))) 

make_matrix = 
  make_df |> 
  filter( !is.na(Code) ) |> 
  select(`Agriculture, forestry, fishing, and hunting`:`Noncomparable imports and rest-of-the-world adjustment /2/`) |> 
  as.matrix()

commodity_output_col = 
  make_df |> 
  filter( Name == "Total Commodity Output") |> 
  select(`Agriculture, forestry, fishing, and hunting`:`Noncomparable imports and rest-of-the-world adjustment /2/`) |> 
  unlist()

q_hat = 
  diag(commodity_output_col)
q_hat.inv =
  solve(q_hat)

industry_output_col = 
  make_df$`Total Industry Output`[!is.na(make_df$Code)]

g_hat = diag(industry_output_col)
g_hat.inv = solve(g_hat)

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
#    "data/raw/input-output-tables/sector-level/com-com_total_requirement.xlsx",
#    skip = 6
#  ) |> 
#  filter( !is.na(Code) ) |> 
#  select(!c(Code,`Commodity Description`)) |> 
#  as.matrix()
#
#error_matrix = 
#  comp_total_req_matrix[-c(400:401),-c(400:401)] - total_req_matrix

ind_ind_total_req_matrix = 
  direct_coef_matrix%*%total_req_matrix
