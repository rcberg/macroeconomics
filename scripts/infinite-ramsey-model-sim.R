
# note: uses CRRA utility function of the form c^(1-g) /(1-g)
# and a cobb-douglas production function of the (per-capita) form f^a

##############################

a = 0.33 # capital share of output
p = 0.02 # discount rate
g = 1 # crra utility parameter
d = 0.015 # depreciation rate
n = 0.01 # workforce/population growth rate
TE = 300 # terminal time period. needs to be big enough for model to converge, small enough to run up to 5k iterations.

##############################

u_c = 
  function(c, g){
    if( g == 1 ){
      u = log(c)
    }else{
      u = (c^(1-g))/(1-g)
    }
    return(u)
  }

u_prime = 
  function(c, g){
    u_prime = c^(-g)
    return(u_prime)
  }

little_f = 
  function(k, a){
    lil_f = k^a
    return(lil_f)
  }

little_f_prime =
  function(k, a){
    lil_f_prime = a*k^(a-1)
    return(lil_f_prime)
  }

k_next = 
  function(k,c, a, n, d){
    
    k_next = k + little_f(k, a) - (n + d)*k - c 
    k_out = max(k_next,0)
    
    return(k_out)
  }

c_next = 
  function(k,c, a, n, d, g, p){
    
    c_next = c + c*(little_f_prime(k,a) - p - d)/g
    c_out = max(c_next,0)
    
    return(c_out)
  }

stationary_k_curve = 
  function(k, a, n, d){
    
    c = k^a - (n+d)*k
    
    return(c)
  }
    
k_ss = 
  (a/(p + d))^(1/(1-a))

c_ss = 
  k_ss^a - (n+d)*k_ss

shooter_function = 
  function(c0, k0, tmax = TE){
    
    if( c0 > little_f(k0,a) ){stop("Consumption is infeasible.")}
    
    c_vec = rep(0, tmax )
    k_vec = rep(0, tmax )
    
    c_vec[1] = c0
    k_vec[1] = k0
    
    for(t in 2:tmax){
      
      
      k_vec[t] = 
        k_next( k_vec[(t-1)], c_vec[(t-1)], a, n, d)
      
      c_vec[t] = 
        c_next( k_vec[(t-1)], c_vec[(t-1)], a, n, d, g, p)
      
    }
    
    return(list(k_vec,c_vec))
  }

bisection_function = 
  function( k0, tmax = TE, tol = 1e-4, max_iter = 2000, k_ter = k_ss, c_ter = c_ss
           ){
    
    c0_upper = little_f(k0, a)
    
    if( k0 > k_ter ){
      
      c0_lower = stationary_k_curve(k0, a, n, d)
    
    }else{
      c0_lower = 0
    }
    
    i=0
    error_k = tol + 1
    error_c = tol + 1
    while(i < max_iter & (error_k > tol | error_c > tol) ){
      
      c0 = (c0_lower+c0_upper)/2
      
      k_vec = shooter_function(c0, k0)[[1]]
      c_vec = shooter_function(c0, k0)[[2]]
      
      if(k_vec[tmax] > k_ter & c_vec[tmax] < c_ter){
        c0_lower = c0 
      }else if(k_vec[tmax] < k_ter & c_vec[tmax] < c_ter){
        c0_upper = c0
      }else if(k_vec[tmax] < k_ter & c_vec[tmax] > c_ter){
        c0_upper = c0
      }
      
      error_k = abs(k_vec[tmax] - k_ter)
      error_c = abs(c_vec[tmax] - c_ter)
      
      i = i+1
      
      }
      
    out_tbl = 
      data.frame(
        capital = k_vec,
        consump = c_vec,
        time = 1:tmax,
        iter = rep(i, TE),
        error_cap = rep(error_k, TE),
        error_cons = rep(error_c, TE)
      )
    
    return(out_tbl)
    }

# heads up; convergence with this shooting algorithm is VERY sensitive to (n, d, p) choices!!
# if you choose different (n, d, p) may not converge before the maximum iterations, and may actually diverge

lower_sp = 
  bisection_function(0.2*k_ss)

upper_sp = 
  bisection_function(2*k_ss)

# plots; let's use ggplot2 to make these look smooth

library(tidyverse)

ggplot( ) + 
  geom_line(
    data = lower_sp |> 
      bind_rows(
        upper_sp
      ),
    aes(x = capital, y = consump),
    linewidth=1,color='red') + 
  geom_vline(xintercept = k_ss, linewidth = 1, color = 'navyblue') +
  geom_line(data = 
              tibble(
                k = seq(0,(n+d)^(-(1/(1-a))),1),
                c = stationary_k_curve(k,a,n,d)
              ),
            aes( x = k, y = c),
            linewidth = 1
            ) +
  labs(title = paste0("Stationary path to steady state"), x = "Capital per worker", y = "Consumption per worker") + 
  theme_minimal()

lower_sp  |> 
  ggplot(aes(x = time, y = consump) ) + 
  geom_line(linewidth = 1) + 
  labs(title = "Consumption", x = "Model time", y = "Output units") + 
  theme_minimal()

lower_sp  |> 
  mutate( utility_cur = u_c(consump, g)) |> 
  ggplot(aes(x = time, y = utility_cur) ) + 
  geom_line(linewidth = 1) + 
  labs(title = "Utility (current)", x = "Model time", y = "Utility units") + 
  theme_minimal()

lower_sp  |> 
  ggplot(aes(x = time, y = capital) ) + 
  geom_line(linewidth = 1) + 
  labs(title = "Capital", x = "Model time", y = "Capital units") + 
  theme_minimal()

lower_sp  |> 
  ggplot(aes(x = time, y = (little_f(capital,a)-consump)/little_f(capital,a) ) ) + 
  geom_line(linewidth = 1) + 
  labs(title = "Savings rate", x = "Model time", y = "Output units") + 
  theme_minimal()

lower_sp  |> 
  ggplot(aes(x = time, y = u_prime(consump,g)) ) + 
  geom_line(linewidth = 1) + 
  labs(title = "Lagrange Multiplier", x = "Model time", y = "Utility units") + 
  theme_minimal()

lower_sp |> 
  ggplot(aes(x = time, y = little_f_prime(capital,a)) ) + 
  geom_line(linewidth = 1) + 
  labs(title = "Rate of return", x = "Model time", y = "Output units") + 
  theme_minimal()
