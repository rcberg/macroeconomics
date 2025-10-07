function f!(F, z, v)
    # see Albouy and Stuart (2020) Supp. Appendix B and C for details
    
    u, w, r, p, x_n, y_n, x, y, n_x, n_y, l_x, l_y, k_x, k_y, k = z.^2 # squaring the inbound argument vector acts as a nonnegativity constraint
    s_y, θ_l, θ_n, ϕ_l, ϕ_n, L_tot, N_tot, ρ, σ, τ, q, a_x, a_y = v
    
    T = τ*( w + r*L_tot/N_tot + ρ*k/N_tot ) # taxes
    F[1] = (u/q)*( (1-s_y)^σ + (s_y^σ)*p^(1-σ) )^(1/(1-σ)) - ((1-τ)*( w + r*L_tot/N_tot + ρ*k/N_tot ) + T) # expenditure minimization
    F[2] = ( (θ_l^σ)*(r^(1-σ)) + (θ_n^σ)*(w^(1-σ)) + ((1-θ_l-θ_n)^σ)*(ρ^(1-σ)) )^(1/(1-σ)) - a_x # traded good p=mc
    F[3] = ( (ϕ_l^σ)*r^(1-σ) + (ϕ_n^σ)*w^(1-σ) + ((1-ϕ_l-ϕ_n)^σ)*(ρ^(1-σ)) )^(1/(1-σ)) - p*a_y # home good p=mc
    F[4] = x_n + p*y_n - ((1-τ)*( w + r*L_tot/N_tot + ρ*k/N_tot ) + T) # walras law
    F[5] = (s_y*(y_n^(-1/σ)))/((1-s_y)*(x_n^(-1/σ))) - p # indifference condition
    F[6] = (θ_n^(σ))*(w^(-σ))*( (θ_l^σ)*r^(1-σ) + (θ_n^σ)*w^(1-σ) + ((1-θ_l-θ_n)^σ)*(ρ^(1-σ)) )^(σ/(1-σ)) - a_x*n_x/x # traded good labor demand
    F[7] = (θ_l^(σ))*(r^(-σ))*( (θ_l^σ)*r^(1-σ) + (θ_n^σ)*w^(1-σ) + ((1-θ_l-θ_n)^σ)*(ρ^(1-σ)) )^(σ/(1-σ)) - a_x*l_x/x # traded good land demand
    F[8] = ((1-θ_l-θ_n)^(σ))*(ρ^(-σ))*( (θ_l^σ)*r^(1-σ) + (θ_n^σ)*w^(1-σ) + ((1-θ_l-θ_n)^σ)*(ρ^(1-σ)) )^(σ/(1-σ)) - a_x*k_x/x # traded good capital demand
    F[9] = (ϕ_n^σ)*(w^(-1/σ))*( (ϕ_l^σ)*r^(1-σ) + (ϕ_n^σ)*w^(1-σ) + ((1-ϕ_l-ϕ_n)^σ)*(ρ^(1-σ)) )^(σ/(1-σ)) - a_y*n_y/y # home good labor demand
    F[10] = (ϕ_l^σ)*(r^(-1/σ))*( (ϕ_l^σ)*r^(1-σ) + (ϕ_n^σ)*w^(1-σ) + ((1-ϕ_l-ϕ_n)^σ)*(ρ^(1-σ)) )^(σ/(1-σ)) - a_y*l_y/y # home good land demand
    F[11] = ((1-ϕ_l-ϕ_n)^σ)*(ρ^(-1/σ))*( (ϕ_l^σ)*r^(1-σ) + (ϕ_n^σ)*w^(1-σ) + ((1-ϕ_l-ϕ_n)^σ)*(ρ^(1-σ)) )^(σ/(1-σ)) - a_y*k_y/y # home good capital demand
    F[12] = N_tot - n_x - n_y # labor market clearing
    F[13] = L_tot - l_x - l_y # land market clearing
    F[14] = k - k_x - k_y # capital market clearing
    F[15] = y - y_n*N_tot # home good market clearing
    return nothing
end 