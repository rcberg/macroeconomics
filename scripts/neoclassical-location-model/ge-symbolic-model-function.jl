using Symbolics

Symbolics.@variables u w r p x y s_y θ_l θ_n ϕ_l ϕ_n ρ σ q

u_xy = q*( (1-s_y)*x^((σ - 1)/σ) + s_y*y^((σ - 1)/σ) )^(σ/(σ - 1))
e_pu = (u/q)*( (1-s_y)^σ + (s_y^σ)*p^(1-σ) )^(1/(1-σ))
c_x_prw = ( (θ_l^σ)*(r^(1-σ)) + (θ_n^σ)*(w^(1-σ)) + ((1-θ_l-θ_n)^σ)*(ρ^(1-σ)) )^(1/(1-σ))
c_y_prw = ( (ϕ_l^σ)*r^(1-σ) + (ϕ_n^σ)*w^(1-σ) + ((1-ϕ_l-ϕ_n)^σ)*(ρ^(1-σ)) )^(1/(1-σ))

utility_function = build_function(u_xy, q, x, y, s_y, σ, expression = false)
expenditure_function = build_function(e_pu, u, q, p, s_y, σ, expression = false )
cost_function_x = build_function(c_x_prw, r, w, ρ, θ_l, θ_n, σ, expression = false)
cost_function_y = build_function(c_y_prw, r, w, ρ, ϕ_l, ϕ_n, σ, expression = false)

D_w = Differential(w)
D_r = Differential(r)
D_i = Differential(ρ)
D_x = Differential(x)
D_y = Differential(y)

dc_x_w = expand_derivatives(D_w(c_x_prw))
dc_y_w = expand_derivatives(D_w(c_y_prw))
dc_x_r = expand_derivatives(D_r(c_x_prw))
dc_y_r = expand_derivatives(D_r(c_y_prw))
dc_x_i = expand_derivatives(D_i(c_x_prw))
dc_y_i = expand_derivatives(D_i(c_y_prw))
du_x = expand_derivatives(D_x(u_xy))
du_y = expand_derivatives(D_y(u_xy))

foc_c_x_w = build_function(dc_x_w, r, w, ρ, θ_l, θ_n, σ, expression = false)
foc_c_x_r = build_function(dc_x_r, r, w, ρ, θ_l, θ_n, σ, expression = false)
foc_c_x_i = build_function(dc_x_i, r, w, ρ, θ_l, θ_n, σ, expression = false)
foc_c_y_w = build_function(dc_y_w, r, w, ρ, ϕ_l, ϕ_n, σ, expression = false)
foc_c_y_r = build_function(dc_y_r, r, w, ρ, ϕ_l, ϕ_n, σ, expression = false)
foc_c_y_i = build_function(dc_y_i, r, w, ρ, ϕ_l, ϕ_n, σ, expression = false)
foc_dx = build_function(du_x, q, x, y, s_y, σ, expression = false)
foc_dy = build_function(du_y, q, x, y, s_y, σ, expression = false)

function f!(F, z, v)
    # large city. see Albouy and Stuart (2020) Supp. Appendix B and C for details
    sig = v.σ 
    
    u, w, r, p, x_n, y_n, x, y, n_x, n_y, l_x, l_y, k_x, k_y, k = z.^2 # squaring the inbound argument vector acts as a nonnegativity constraint
    s_y, θ_l, θ_n, ϕ_l, ϕ_n, L_tot, N_tot, ρ, σ, τ, q, a_x, a_y = v
    
    T = τ*( w + r*L_tot/N_tot + ρ*k/N_tot ) # taxes
    F[1] = expenditure_function(u, q, p, s_y, σ) - ((1-τ)*( w + r*L_tot/N_tot + ρ*k/N_tot ) + T) # expenditure minimization
    F[2] = cost_function_x(r, w, ρ, θ_l, θ_n, σ) - a_x # traded good p=mc
    F[3] = cost_function_y(r, w, ρ, ϕ_l, ϕ_n, σ) - p*a_y # home good p=mc
    F[4] = x_n + p*y_n - ((1-τ)*( w + r*L_tot/N_tot + ρ*k/N_tot ) + T) # walras law
    F[5] = foc_dy(q, x, y, s_y, σ)/foc_dx(q, x, y, s_y, σ) - p # indifference condition
    F[6] = foc_c_x_w(r, w, ρ, θ_l, θ_n, σ) - a_x*n_x/x # traded good labor demand
    F[7] = foc_c_x_r(r, w, ρ, θ_l, θ_n, σ) - a_x*l_x/x # traded good land demand
    F[8] = foc_c_x_i(r, w, ρ, θ_l, θ_n, σ) - a_x*k_x/x # traded good capital demand
    F[9] = foc_c_y_w( r, w, ρ, ϕ_l, ϕ_n, σ) - a_y*n_y/y # home good labor demand
    F[10] = foc_c_y_r( r, w, ρ, ϕ_l, ϕ_n, σ) - a_y*l_y/y # home good land demand
    F[11] = foc_c_y_i( r, w, ρ, ϕ_l, ϕ_n, σ) - a_y*k_y/y # home good capital demand
    F[12] = N_tot - n_x - n_y # labor market clearing
    F[13] = L_tot - l_x - l_y # land market clearing
    F[14] = k - k_x - k_y # capital market clearing
    F[15] = y - y_n*N_tot # home good market clearing
    return nothing
end 

function g!(G, z, v)
    # small city
    
    w, r, p, x_n, y_n, x, y, n_x, n_y, l_x, l_y, k_x, k_y, N_tot, k = z.^2 # squaring the inbound argument vector acts as a nonnegativity constraint
    s_y, θ_l, θ_n, ϕ_l, ϕ_n, L_tot, u, ρ, σ, τ, q, a_x, a_y, Ltot_rw, Ntot_rw, k_nw = v
    
    T = τ*( w + r*Ltot_rw/Ntot_rw + ρ*k_nw/Ntot_rw ) # taxes
    G[1] = expenditure_function(u, q, p, s_y, σ) - ((1-τ)*( w + r*Ltot_rw/Ntot_rw + ρ*k_nw/Ntot_rw ) + T) # expenditure minimization
    G[2] = cost_function_x(r, w, ρ, θ_l, θ_n, σ) - a_x # traded good p=mc
    G[3] = cost_function_y(r, w, ρ, ϕ_l, ϕ_n, σ) - p*a_y # home good p=mc
    G[4] = x_n + p*y_n - ((1-τ)*( w + r*L_tot/N_tot + ρ*k/N_tot ) + T) # walras law
    G[5] = foc_dy(q, x, y, s_y, σ)/foc_dx(q, x, y, s_y, σ) - p # indifference condition
    G[6] = foc_c_x_w(r, w, ρ, θ_l, θ_n, σ) - a_x*n_x/x # traded good labor demand
    G[7] = foc_c_x_r(r, w, ρ, θ_l, θ_n, σ) - a_x*l_x/x # traded good land demand
    G[8] = foc_c_x_i(r, w, ρ, θ_l, θ_n, σ) - a_x*k_x/x # traded good capital demand
    G[9] = foc_c_y_w( r, w, ρ, ϕ_l, ϕ_n, σ) - a_y*n_y/y # home good labor demand
    G[10] = foc_c_y_r( r, w, ρ, ϕ_l, ϕ_n, σ) - a_y*l_y/y # home good land demand
    G[11] = foc_c_y_i( r, w, ρ, ϕ_l, ϕ_n, σ) - a_y*k_y/y # home good capital demand
    G[12] = N_tot - n_x - n_y # labor market clearing
    G[13] = L_tot - l_x - l_y # land market clearing
    G[14] = k - k_x - k_y # capital market clearing
    G[15] = y - y_n*N_tot # home good market clearing
    return nothing
end