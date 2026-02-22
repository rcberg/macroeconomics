# NOTE: needs Symbolics.jl, NonlinearSolve.jl, and LinearAlgebra.jl
# normally just import those at the top of whatever project uses these functions...
# ...but b/c 'gradient()' demands a qualified name, we must import Symbolics.jl or the function is DoA
using Symbolics

# the following 2 functions bound our variables >= 0 via transformation of variables.
# `sigmoid()` constrains the solution to be >=0 and <= [limit]. `logit()` is just the inverse.
function sigmoid(val, limit)
     limit / (1.0 + exp(-val))
end

function logit(val, limit)
    if val isa Number
        ϵ = 1e-12
        x_safe = clamp(val, ϵ, limit - ϵ) # `clamp` is interesting, ?clamp in REPL to learn more.
        return log(x_safe / (limit - x_safe) )
    else
        return log(val / (limit - val))
    end
end

function economy_model( parameters, vars_guess )
    p_vals = [ 
        parameters.γ_L,
        parameters.γ_N,
        parameters.ρ_L,
        parameters.ρ_N,
        parameters.η_X,
        parameters.η_Y,
        parameters.q,
        parameters.a_x,
        parameters.a_y,
        parameters.L,
        parameters.σ_D,
        parameters.σ_X,
        parameters.σ_Y,
        parameters.N_T,
        parameters.ι]

    guessvec = [
        vars_guess.x,
        vars_guess.y,
        vars_guess.p,
        vars_guess.w,
        vars_guess.r,
        vars_guess.λ,
        vars_guess.L_x,
        vars_guess.L_y,
        vars_guess.N,
        vars_guess.N_x,
        vars_guess.N_y,
        vars_guess.K,
        vars_guess.K_x,
        vars_guess.K_y]

    @variables v[1:14]
    @variables q a_x a_y L σ_D σ_X σ_Y N_T ι
    @variables γ_L γ_N ρ_L ρ_N η_X η_Y η_N
    
    # makes it Inf times easier to set up the model using variable and parameter names.
    idx_u = ( 
        x = 1,
        y = 2,
        p = 3,
        w = 4,
        r = 5,
        λ = 6,
        L_x = 7,
        L_y = 8,
        N = 9,
        N_x = 10,
        N_y = 11,
        K = 12,
        K_x = 13,
        K_y = 14
    )
    
    idx_p = (
        γ_L = 1, 
        γ_N = 2, 
        ρ_L = 3, 
        ρ_N = 4, 
        η_X = 5, 
        η_Y = 6, 
        q = 7, 
        a_x = 8, 
        a_y = 9, 
        L = 10, 
        σ_D = 11, 
        σ_X = 12, 
        σ_Y = 13, 
        N_T = 14,
        ι = 15
    )
    
    # this is transformation of variables ffor the symbolic variables/expressions ONLY
    d = (
        x = exp(v[idx_u.x]),
        y = exp(v[idx_u.y]),
        p = exp(v[idx_u.p]),
        w  = exp(v[idx_u.w]),
        r  = exp(v[idx_u.r]),
        λ = exp(v[idx_u.λ]),
        L_x = sigmoid(v[idx_u.L_x], L),
        L_y = sigmoid(v[idx_u.L_y], L),
        N = exp(v[idx_u.N]),
        N_x = exp(v[idx_u.N_x]),
        N_y = exp(v[idx_u.N_y]),
        K = exp(v[idx_u.K]),
        K_x = exp(v[idx_u.K_x]),
        K_y = exp(v[idx_u.K_y])
        )
    
    α = (σ_D - 1)/σ_D
    β = (σ_X - 1)/σ_X
    χ = (σ_Y - 1)/σ_Y
    
    # utility and production/profit functions
    # want to provide the option for cobb-douglas
    if p_vals[idx_p.σ_D] == 1.0
        u_xy = q*(d.x^η_X)*(d.y^η_Y)*((N_T - d.N)^(1 - η_X - η_Y))
    else
        u_xy = q*( η_X*(d.x)^α + η_Y*(d.y)^α + (1 - η_X - η_Y)*(N_T - d.N)^α )^(1/α)
    end

    if p_vals[idx_p.σ_X] == 1.0
        pf_x = a_x*(d.L_x^γ_L)*(d.N_x^γ_N)*(d.K_x^(1 - γ_L - γ_N))
    else
        pf_x = a_x*( γ_L*(d.L_x)^β + γ_N*(d.N_x)^β + (1 - γ_L - γ_N)*(d.K_x)^β )^(1/β)
    end
    
    if p_vals[idx_p.σ_Y] == 1.0
        pf_y = a_y*(d.L_y^ρ_L)*(d.N_y^ρ_N)*(d.K_y^(1 - ρ_L - ρ_N)) 
    else
        pf_y = a_y*( ρ_L*(d.L_y)^χ + ρ_N*(d.N_y)^χ + (1 - ρ_L - ρ_N)*(d.K_y)^χ )^(1/χ)
    end

    Λ_u = u_xy + d.λ * ( d.w*d.N + d.r*L + ι*d.K - d.x - d.p*d.y )
    Π_x = pf_x - ( d.w*d.N_x + d.r*d.L_x + ι*d.K_x )
    Π_y = (d.p)*pf_y - ( d.w*d.N_y + d.r*d.L_y + ι*d.K_y )
    
    # set up the system of equations with FOCs and market clearing conditions
    h_focs = Symbolics.gradient(Λ_u, [v[idx_u.x], v[idx_u.y], v[idx_u.N], v[idx_u.λ]] )
    firm_x_focs = Symbolics.gradient(Π_x, [v[idx_u.L_x], v[idx_u.N_x], v[idx_u.K_x]] )
    firm_y_focs = Symbolics.gradient(Π_y, [v[idx_u.L_y], v[idx_u.N_y], v[idx_u.K_y]] )
    
    #k_clearing = d.K - d.K_x - d.K_y # fixed ι means dropping 1 equation from the system. K worked best.
    l_clearing = L - d.L_x - d.L_y
    n_clearing = d.N - d.N_x - d.N_y
    
    x_clearing = d.x - pf_x
    y_clearing = d.y - pf_y
    
    system_eqs = [firm_x_focs; firm_y_focs; h_focs; y_clearing; x_clearing; l_clearing; n_clearing]
    
    params = [γ_L, γ_N, ρ_L, ρ_N, η_X, η_Y, q, a_x, a_y, L, σ_D, σ_X, σ_Y, N_T, ι]
    full_sys = Symbolics.expand(system_eqs)
    
    # we supply out own Jacobian to make it easier for the solver
    j_matrix = Symbolics.jacobian(full_sys, v)
    
    # these 2 lines turn our symbolic functions into "regular" julia functions, for `solve()`
    f_ge! = build_function(full_sys, v, params, expression=false)[2]
    j_ip! = build_function(j_matrix, v, params, expression=false)[2]
    
    ge_system = NonlinearFunction(f_ge!; jac = j_ip!)

    # transforming the actual, supplied guesses and feeding it to the solver...
    u0 = [ log.(guessvec[1:6]); logit.(guessvec[7:8], p_vals[idx_p.L]); log.(guessvec[9:14]) ]
    problem = NonlinearProblem( ge_system, u0, p_vals )
    
    sol = solve(problem, LevenbergMarquardt(), abstol = 1e-9, reltol = 1e-12)
    
    raw_solution = sol.u
    if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.StalledSuccess
        successful = 1
    else
        successful = 0
    end
    # ...and transforming the output back into correct, zero-bounded zolutions to the model
    solution = [ exp.(raw_solution[1:6]); sigmoid.(raw_solution[7:8], p_vals[idx_p.L]); exp.(raw_solution[9:14]) ]
    
    alfa = (p_vals[idx_p.σ_D] - 1.0)/p_vals[idx_p.σ_D]
    beta = (p_vals[idx_p.σ_X] - 1.0)/p_vals[idx_p.σ_X]
    chi =  (p_vals[idx_p.σ_Y] - 1.0)/p_vals[idx_p.σ_Y]
    
    profit_x = solution[idx_u.x] - ( solution[idx_u.w]*solution[idx_u.N_x] + solution[idx_u.r]*solution[idx_u.L_x] + p_vals[idx_p.ι]*solution[idx_u.K_x] )
    profit_y = solution[idx_u.p]*solution[idx_u.y] - ( solution[idx_u.w]*solution[idx_u.N_y] + solution[idx_u.r]*solution[idx_u.L_y] + p_vals[idx_p.ι]*solution[idx_u.K_y] )
    profit_total = profit_x + profit_y
    
    # gradient check
    # these should all be tiny on the order of << 1e-10
    walras_law = (solution[idx_u.w]*solution[idx_u.N] + solution[idx_u.r]*p_vals[idx_p.L] + p_vals[idx_p.ι]*solution[idx_u.K] - solution[idx_u.x] - solution[idx_u.p]*solution[idx_u.y])
    l_mkt_clearing = (p_vals[idx_p.L] - solution[idx_u.L_x] - solution[idx_u.L_y]) # land market clearing
    n_mkt_clearing = (solution[idx_u.N] - solution[idx_u.N_x] - solution[idx_u.N_y]) # labor market clearing 
    k_mkt_clearing = (solution[idx_u.K] - solution[idx_u.K_x] - solution[idx_u.K_y]) # capital market clearing
    good_x_profit = profit_x # firm X zero-profit condition
    good_y_profit = profit_y # firm Y zero-profit condition
    industrywide_profit = profit_total # total zero profit condition
    λ_multiplier = solution[idx_u.λ] - p_vals[idx_p.η_X]  * (solution[idx_u.x]^(  alfa - 1)) * p_vals[idx_p.q]   * ( p_vals[idx_p.η_X]*(solution[idx_u.x])^alfa   + p_vals[idx_p.η_Y]*(solution[idx_u.y])^alfa   + (1 - p_vals[idx_p.η_X] - p_vals[idx_p.η_Y])*(p_vals[idx_p.N_T] - solution[idx_u.N])^alfa )^((1/alfa) - 1)  # consumption foc
    good_x_w_foc = solution[idx_u.w] - p_vals[idx_p.γ_N]  * (solution[idx_u.N_x]^(beta - 1)) * p_vals[idx_p.a_x] * ( p_vals[idx_p.γ_L]*(solution[idx_u.L_x])^beta + p_vals[idx_p.γ_N]*(solution[idx_u.N_x])^beta + (1 - p_vals[idx_p.γ_L] - p_vals[idx_p.γ_N])*(solution[idx_u.K_x])^beta )^((1/beta) - 1) # good x wage foc
    good_y_w_foc = solution[idx_u.w] - solution[idx_u.p] * p_vals[idx_p.ρ_N]  * (solution[idx_u.N_y]^(chi  - 1)) * p_vals[idx_p.a_y] * ( p_vals[idx_p.ρ_L]*(solution[idx_u.L_y])^chi  + p_vals[idx_p.ρ_N]*(solution[idx_u.N_y])^chi  + (1 - p_vals[idx_p.ρ_L] - p_vals[idx_p.ρ_N])*(solution[idx_u.K_y])^chi  )^((1/chi) -  1) # good y wage foc
    good_x_r_foc = solution[idx_u.r] - p_vals[idx_p.γ_L]  * (solution[idx_u.L_x]^(beta - 1)) * p_vals[idx_p.a_x] * ( p_vals[idx_p.γ_L]*(solution[idx_u.L_x])^beta + p_vals[idx_p.γ_N]*(solution[idx_u.N_x])^beta + (1 - p_vals[idx_p.γ_L] - p_vals[idx_p.γ_N])*(solution[idx_u.K_x])^beta )^((1/beta) - 1) # good x land foc
    good_y_r_foc = solution[idx_u.r] - solution[idx_u.p] * p_vals[idx_p.ρ_L]  * (solution[idx_u.L_y]^(chi  - 1)) * p_vals[idx_p.a_y] * ( p_vals[idx_p.ρ_L]*(solution[idx_u.L_y])^chi  + p_vals[idx_p.ρ_N]*(solution[idx_u.N_y])^chi  + (1 - p_vals[idx_p.ρ_L] - p_vals[idx_p.ρ_N])*(solution[idx_u.K_y])^chi  )^((1/chi) -  1) # good y land foc
    good_x_k_foc = p_vals[idx_p.ι]   - (1 - p_vals[idx_p.γ_L] - p_vals[idx_p.γ_N]) * (solution[idx_u.K_x]^(beta - 1)) * p_vals[idx_p.a_x] * ( p_vals[idx_p.γ_L]*(solution[idx_u.L_x])^beta + p_vals[idx_p.γ_N]*(solution[idx_u.N_x])^beta + (1 - p_vals[idx_p.γ_L] - p_vals[idx_p.γ_N])*(solution[idx_u.K_x])^beta )^((1/beta) - 1) # good x capt foc
    good_y_k_foc = p_vals[idx_p.ι]   - solution[idx_u.p] * (1 - p_vals[idx_p.ρ_L] - p_vals[idx_p.ρ_N]) * (solution[idx_u.K_y]^(chi  - 1)) * p_vals[idx_p.a_y] * ( p_vals[idx_p.ρ_L]*(solution[idx_u.L_y])^chi  + p_vals[idx_p.ρ_N]*(solution[idx_u.N_y])^chi  + (1 - p_vals[idx_p.ρ_L] - p_vals[idx_p.ρ_N])*(solution[idx_u.K_y])^chi  )^((1/chi) -  1) # good y capt foc
    
    con_vec = [walras_law; l_mkt_clearing; n_mkt_clearing; k_mkt_clearing; good_x_profit; good_y_profit; industrywide_profit; λ_multiplier; good_x_w_foc; good_y_w_foc; good_x_r_foc; good_y_r_foc; good_x_k_foc; good_y_k_foc]
    sol_vec = solution

    return(con_vec, sol_vec, successful)
end