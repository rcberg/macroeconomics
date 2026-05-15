using JuMP, Ipopt

function fill_lower_tri!(H, func, args, n)
        k = 1
        for col in 1:n
            for row in col:n
                H[row, col] = func[k](args)
                k += 1
            end
        end
    end

function build_model(fns; 
    Q, A_X, A_Y, L_0, 
    u_bar = nothing, N_fixed = nothing, 
    ηx_fixed = nothing, γL_fixed = nothing, γN_fixed = nothing, ρL_fixed = nothing, ρN_fixed = nothing,
    R_fixed = nothing, I_fixed = nothing, T_fixed = nothing,
    τ, ι_bar,
    λ_L = 0.17, λ_N = 0.7)

    ϵ = 0.0

    # parameter targets from albouy and stuart's table 1
    target_sy = 0.36
    target_θL = 0.025 
    target_θN = 0.825
    target_ϕL = 0.233
    target_ϕN = 0.617

    # initialize the optimizer
    urban_economy = Model(Ipopt.Optimizer)
    #set_optimizer_attribute(urban_economy, "print_level", 0) # un-comment to silence output
    
    # incorporate the compiled cost/expenditure functions, including gradients, jacobians, and hessians
    # note how the jacobian and hessian are supplied for the optimization algorithm.
    # traded good (x)
    @operator(urban_economy, op_cX, 4,
    (r, w, γ_L, γ_N) -> fns.cX_fun([r, w, γ_L, γ_N]),
    (g, r, w, γ_L, γ_N) -> fns.dcX_fun_ip!(g, [r, w, γ_L, γ_N]),
    (H, r, w, γ_L, γ_N) -> fill_lower_tri!(H, fns.d2cX_hess_fns, [r, w, γ_L, γ_N], 4)
    )
    @operator(urban_economy, op_∂cX_∂r, 4,
    (r, w, γ_L, γ_N) -> fns.dcX_funs[1]([r, w, γ_L, γ_N]),
    (g, r, w, γ_L, γ_N) -> fns.d2cX_row_ip![1](g, [r, w, γ_L, γ_N]),
    (H, r, w, γ_L, γ_N) ->  fill_lower_tri!(H, fns.d3cX_r_hess_fns, [r, w, γ_L, γ_N], 4)
    )
    @operator(urban_economy, op_∂cX_∂w, 4,
    (r, w, γ_L, γ_N) -> fns.dcX_funs[2]([r, w, γ_L, γ_N]),
    (g, r, w, γ_L, γ_N) -> fns.d2cX_row_ip![2](g, [r, w, γ_L, γ_N]),
    (H, r, w, γ_L, γ_N) ->  fill_lower_tri!(H, fns.d3cX_w_hess_fns, [r, w, γ_L, γ_N], 4)
    )
    @operator(urban_economy, op_∂cX_∂i, 4,
    (r, w, γ_L, γ_N) -> fns.dcX_di_fun([r, w, γ_L, γ_N]),
    (g, r, w, γ_L, γ_N) -> fns.d2cX_di_ip!(g, [r, w, γ_L, γ_N]),
    (H, r, w, γ_L, γ_N) -> fill_lower_tri!(H, fns.d3cX_di_hess_fns, [r, w, γ_L, γ_N], 4)
    )
    # home good (y)
    @operator(urban_economy, op_cY, 4,
    (r, w, ρ_L, ρ_N) -> fns.cY_fun([r, w, ρ_L, ρ_N]),
    (g, r, w, ρ_L, ρ_N) -> fns.dcY_fun_ip!(g, [r, w, ρ_L, ρ_N]),
    (H, r, w, ρ_L, ρ_N) -> fill_lower_tri!(H, fns.d2cY_hess_fns, [r, w, ρ_L, ρ_N], 4)
    )
    @operator(urban_economy, op_∂cY_∂r, 4,
    (r, w, ρ_L, ρ_N) -> fns.dcY_funs[1]([r, w, ρ_L, ρ_N]),
    (g, r, w, ρ_L, ρ_N) -> fns.d2cY_row_ip![1](g, [r, w, ρ_L, ρ_N]),
    (H, r, w, ρ_L, ρ_N) -> fill_lower_tri!(H, fns.d3cY_r_hess_fns, [r, w, ρ_L, ρ_N], 4)
    )
    @operator(urban_economy, op_∂cY_∂w, 4,
    (r, w, ρ_L, ρ_N) -> fns.dcY_funs[2]([r, w, ρ_L, ρ_N]),
    (g, r, w, ρ_L, ρ_N) -> fns.d2cY_row_ip![2](g, [r, w, ρ_L, ρ_N]),
    (H, r, w, ρ_L, ρ_N) -> fill_lower_tri!(H, fns.d3cY_w_hess_fns, [r, w, ρ_L, ρ_N], 4)
     )
    @operator(urban_economy, op_∂cY_∂i, 4,
    (r, w, ρ_L, ρ_N) -> fns.dcY_di_fun([r, w, ρ_L, ρ_N]),
    (g, r, w, ρ_L, ρ_N) -> fns.d2cY_di_ip!(g, [r, w, ρ_L, ρ_N]),
    (H, r, w, ρ_L, ρ_N) -> fill_lower_tri!(H, fns.d3cY_di_hess_fns, [r, w, ρ_L, ρ_N], 4)
     )
    # expenditure
    @operator(urban_economy, op_e, 2,
    (p, η_x) -> fns.exp_fun([p, η_x]),
    (g, p, η_x) -> fns.dexp_fun_ip!(g,[ p, η_x]),
    (H, p, η_x) -> begin
        args = [p, η_x]
        H[1, 1] = fns.d2e_hess_fns[1](args)
        H[2, 1] = fns.d2e_hess_fns[2](args)
        H[2, 2] = fns.d2e_hess_fns[3](args)
        end 
    )
    @operator(urban_economy, op_∂e_∂p, 2,
    (p, η_x) -> fns.dexp_dp_fn([p, η_x]),
    (g, p, η_x) -> fns.d2exp_dp_ip!(g, [p, η_x]) ,
    (H, p, η_x) -> begin
        args = [p, η_x]
        H[1, 1] = fns.d3e_dp_hess_fns[1](args)
        H[2, 1] = fns.d3e_dp_hess_fns[2](args)
        H[2, 2] = fns.d3e_dp_hess_fns[3](args)
        end 
    )

    # this basically checks whether to solve the large or small city model
    # necessary because small city takes certain large city values as exogenous
    if isnothing(u_bar) == isnothing(N_fixed)

        error("Make sure u_bar and N_fixed are either both nothing or both Floats.")

    elseif isnothing(u_bar)

        L_total = L_0
        @variable(urban_economy, u_bar_v, start = 1.0)
        @variable(urban_economy, N_total == N_fixed)
        @variable(urban_economy, NX >= ϵ, start = 0.5*N_fixed)
        @variable(urban_economy, NY >= ϵ, start = 0.5*N_fixed)
        @variable(urban_economy, 0.0 <= η_x <= 1.0, start = 0.5 )
        @variable(urban_economy, 0.0 <= γ_L <= 1.0, start = 0.3 )
        @variable(urban_economy, 0.0 <= γ_N <= 1.0, start = 0.3 )
        @variable(urban_economy, 0.0 <= ρ_L <= 1.0, start = 0.3 )
        @variable(urban_economy, 0.0 <= ρ_N <= 1.0, start = 0.3 )
        @variable(urban_economy, R >= ϵ, start = 0.1)
        @variable(urban_economy, I >= ϵ, start = 0.1)
        @variable(urban_economy, T >= 0.0, start = 0.1) # needs to be 0.0 to accommodate τ = 0 option
    
        @constraint(urban_economy, γ_L + γ_N <= 1.0 )
        @constraint(urban_economy, ρ_L + ρ_N <= 1.0 )

    else
        if any(isnothing, (ηx_fixed, γL_fixed, γN_fixed, ρL_fixed, ρN_fixed,
            R_fixed, I_fixed, T_fixed))
            error("Small city requires calibrated share parameters and income/transfers")
        end

        L_total = L_0*1e-6
        @variable( urban_economy, u_bar_v == u_bar)
        @variable( urban_economy, N_total >= 0.0, start = L_total)
        @variable(urban_economy, NX >= ϵ, start = 0.5*L_total)
        @variable(urban_economy, NY >= ϵ, start = 0.5*L_total)
        η_x = ηx_fixed
        γ_L = γL_fixed
        γ_N = γN_fixed
        ρ_L = ρL_fixed
        ρ_N = ρN_fixed
        R = R_fixed
        I = I_fixed
        T = T_fixed

    end

    # set up equilibrium quantities
    @variable(urban_economy, w >= ϵ, start = 1.0)
    @variable(urban_economy, p >= ϵ, start = 1.0)
    @variable(urban_economy, r >= ϵ, start = 1.0)
    @variable(urban_economy, x >= ϵ, start = 1.0)
    @variable(urban_economy, y >= ϵ, start = 1.0)
    @variable(urban_economy, X >= ϵ, start = 1.0)
    @variable(urban_economy, Y >= ϵ, start = 1.0)
    @variable(urban_economy, LX >= ϵ, start = 0.5*L_total)
    @variable(urban_economy, LY >= ϵ, start = 0.5*L_total)
    @variable(urban_economy, KX >= ϵ, start = 0.5) 
    @variable(urban_economy, KY >= ϵ, start = 0.5)
    @variable(urban_economy, K >= ϵ, start = 1.0)
    
    if isnothing(u_bar)
        # consumer problem equilibrium equations
        @constraint(urban_economy, R == (r*L_total) / N_total )
        @constraint(urban_economy, I == (ι_bar*K) / N_total )
        @constraint(urban_economy, T == τ*(w + R + I) )
    end

    income = (1 - τ)*(w + R + I) + T

    @constraint(urban_economy, op_e(p, η_x) * u_bar_v / Q == income)
    
    @constraint(urban_economy, op_cX(r, w, γ_L, γ_N) / A_X == 1.0) 
    @constraint(urban_economy, op_cY(r, w, ρ_L, ρ_N) / A_Y == p)
    
    # consumption demands (2a*, 2b*)
    @constraint(urban_economy, x + p * y == income)
    @constraint(urban_economy, op_∂e_∂p(p, η_x) * u_bar_v / Q == y )
    
    # factor demands (3a*-3f*)
    @constraint(urban_economy, op_∂cX_∂w(r, w, γ_L, γ_N) * (X / A_X) == NX )
    @constraint(urban_economy, op_∂cX_∂r(r, w, γ_L, γ_N) * (X / A_X) == LX )
    @constraint(urban_economy, op_∂cX_∂i(r, w, γ_L, γ_N) * (X / A_X) == KX )
    @constraint(urban_economy, op_∂cY_∂w(r, w, ρ_L, ρ_N) * (Y / A_Y) == NY )
    @constraint(urban_economy, op_∂cY_∂r(r, w, ρ_L, ρ_N) * (Y / A_Y) == LY )
    @constraint(urban_economy, op_∂cY_∂i(r, w, ρ_L, ρ_N) * (Y / A_Y) == KY )
    
    # resource constraints & market clearing (4a*-4c*, 6*)
    @constraint(urban_economy, N_total == NX + NY)
    @constraint(urban_economy, L_total == LX + LY)
    @constraint(urban_economy, K == KX + KY)
    @constraint(urban_economy, Y == N_total*y)
    
    # calibrating the large city model
    # objective: minimize distance to table 1 targets
    # uses the square of the sum of the percentage difference
    # (using just the square differences was failing to match the λs.)

    if isnothing(u_bar)
        @objective(urban_economy, Min, 
        ( ((p * y) / (x + p * y) - target_sy)/target_sy )^2 +
        ( ((w * NX) / (1.0 * X) - target_θN)/target_θN )^2 +
        ( ((r * LX) / (1.0 * X) - target_θL)/target_θL )^2 +
        ( ((w * NY) / (p * Y)   - target_ϕN)/target_ϕN )^2 +
        ( ((r * LY) / (p * Y)   - target_ϕL)/target_ϕL )^2 +
        ( (NX/N_total - λ_N) / λ_N)^2 + 
        ( (LX/L_total - λ_L) / λ_L)^2 
        )
    else
        @objective(urban_economy, Min, 0)
    end

    return urban_economy
end
