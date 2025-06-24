
############################### 

max_iter = 1000 # nr of model iterations
a = 0.3 # capital share of income
b = 0.95 # consumer intertemporal discount rate
d = 0.05 # depreciation over the generation cycle
n = 0.01 # population growth rate
g = 0.03 # technology/tfp growth rate
z = n+g+n*g # (1+n)(1+g)
l0 = 0.5 # initial labor
k0 = 0.5 # initial capital
tfp0=0.5 # initial technology/tfp

###############################

# labor and tfp paths

l = l0*(1+n)^(0:max_iter)     # cool R trick a^(vector) is a vector
tfp = tfp0*(1+g)^(0:max_iter)

# steady state and golden rule capital levels
# note production function is cobb-douglas crs in (K,L)

k_ss = ((b*(1-a))/((1+b)*(z+d)))^(1/(1-a))
y_ss = k_ss^a
k_gr = (a/(z+d))^(1/(1-a))

k = c(k0, rep(0,max_iter) )   # initialized capital vector

# capital stock path based on equation of motion

for(t in 1:max_iter){
  
  k[(t+1)] = (k[t]^a)*(b*(1-a))/((1+b)*(1+z)) + k[t]*(1-d)/(1+z)
}

# initialize tax sim

k_tax = 
  c(k_ss, rep(0,max_iter) )

# fun trick i discovered.
# you can get this next result by assuming existence of a golden rule steady state level of k with gov>0, and solve for gov.

gov = (1-a)*k_gr^a - k_gr*((1+b)*(z+d))/b

for(t in 1:max_iter){
  
  k_tax[(t+1)] = (b/((1+b)*(1+z)))*((1-a)*k_tax[t]^a - gov) + k_tax[t]*(1-d)/(1+z)
  
}

library(tidyverse)

baseline_df = 
  tibble(
    time = 0:max_iter,
    capital = k,
    cap_ss = k_ss,
    output = k^a,
    tfp = tfp,
    y_ss = k_ss^a,
    wage = (1-a)*capital^a,
    rate = a*capital^(a-1),
    saving = (b/(1+b))*((1-a)*capital^a)
  ) 

forward_df = 
  tibble(
    time = 0:(2*max_iter+1),
    capital = c(k,k_tax),
    output = (capital)^a,
    wage = (1-a)*capital^a,
    rate = a*capital^(a-1),
    saving = (b/(1+b))*((1-a)*capital^a - gov)
  ) 

## plots ##

baseline_df |> 
  ggplot() + 
  geom_point(
    aes(
      x = time, 
      y = capital
      ),
      color = 'blue'
  ) +
  geom_point(
    aes(
      x = time, 
      y = output
      ),
      color = 'red'
  ) + 
  geom_hline(
    yintercept = k_ss, linetype = 2, color = 'blue', linewidth = 1
  ) + 
  geom_hline(
    yintercept = y_ss, linetype = 2, color = 'red', linewidth = 1
  ) +
  geom_hline(
    yintercept = k_gr, linetype = 2, color = 'green', linewidth = 1
  ) +
  labs( x = "Time", y = "Capital (blue) / Output (red)",
        caption = 'Green line shows Golden Rule level of capital per eff. worker'
        ) + 
  theme_minimal()

# this is cool. you can see how we use the tax to achieve dynamic efficiency

forward_df |> 
  ggplot() + 
  geom_point(
    aes(
      x = time, 
      y = capital
    ),
    color = 'blue'
  ) +
  geom_point(
    aes(
      x = time, 
      y = output
    ),
    color = 'red'
  ) + 
  geom_hline(
    yintercept = k_gr, linetype = 2, color = 'green', linewidth = 1
  ) +
  labs( x = "Time", y = "Capital (blue) / Output (red)"
  ) + 
  theme_minimal()

forward_df |> 
  ggplot() + 
  geom_point(
    aes(
      x = time, 
      y = wage
    ),
    color = 'blue'
  ) +
  labs( x = "Time", y = "Wage (output units)") +
  theme_minimal()

forward_df |> 
  ggplot() + 
  geom_point(
    aes(
      x = time, 
      y = rate
    ),
    color = 'red'
  ) +
  labs( x = "Time", y = "Rate of return") +
  theme_minimal()

# diamond paper diagram 3

cols = c('phi' = 'blue', 'psi' = 'red')
baseline_df |> 
  mutate( phi = a*((b/(1+b))*lag(wage)/(1+z))^(a-1) ,
          psi = a*((b/(1+b))*wage/(1+z))^(a-1) ) |>  
  ggplot(  ) + 
  geom_line(aes( x = phi, y = wage , color = 'phi'), linewidth = 0.75) + 
  geom_line(aes( x = psi, y = wage , color = 'psi'), linewidth = 0.75) +
  scale_color_manual( name = 'Frontier', values = cols) +
  labs(x = 'rate', title = "Diamond (1965) Diagram 3") + 
  theme_minimal()
