using JuMP, Ipopt
using Plots

# 1. Model constants/calibration targets
#σ_D = 0.667
#σ_X = 0.667
#σ_Y = 0.667
#i_bar = 1.0
#L_large = 1000.0
#L_small = L_large*1e-6
#N_large = 1000.0
#τ = 0.0
#λ_L = 0.17
#λ_N = 0.7
# Exogenous Parameters (Fixed for Large City)
#AX = 1.0
#AY = 1.0
#Q = 1.0

function city_function(;AX = 1.0, AY = 1.0, Q = 1.0,
                    N_fixed = nothing, U = nothing,
                    ηx_fixed = nothing, γL_fixed = nothing, γN_fixed = nothing, ρL_fixed = nothing, ρN_fixed = nothing,
                    R_fixed = nothing, I_fixed = nothing, T_fixed = nothing,
                    L_in = 1000.0, τ = 0.0, λ_L = 0.17, λ_N = 0.7 )

    σ_D = 0.667
    σ_X = 0.667
    σ_Y = 0.667
    i_bar = 1.0
    # 2. Initialize large city model
    
    # Using `Ipopt` solver, but others (e.g., KNITRO) could be used if you're licensed.
    urban_city = Model(Ipopt.Optimizer)
    #set_optimizer_attribute(urban_city, "print_level", 0) # Un-comment to silence output
    
    if isnothing(N_fixed) == isnothing(U)
        error("Make sure you supply either U or N.")
    elseif isnothing(U)
        
        N = N_fixed
        L = L_in
        # Calibration Parameters (Variables for Large City)
        @variables(urban_city, begin
        u_bar, (start = 1.0)
        NX, (start = 0.5*N)
        NY, (start = 0.5*N) # Labor Demands
        0 <= η_x <= 1, (start = 0.5)
        0 <= γ_L <= 1, (start = 0.3)
        0 <= γ_N <= 1, (start = 0.3)
        0 <= ρ_L <= 1, (start = 0.3)
        0 <= ρ_N <= 1, (start = 0.3)
        end)
        # nonlabor income and transfers (also large city variables)
        @variables(urban_city, begin
        0 <= R, (start = 0.1)
        0 <= I, (start = 0.1)
        0 <= T, (start = 0.1)
        end)
    else
        if any(isnothing, (ηx_fixed, γL_fixed, γN_fixed, ρL_fixed, ρN_fixed, 
            R_fixed, I_fixed, T_fixed))
            error("Small city requires calibrated share parameters and per-capita income.")
        end

        u_bar = U
        L = L_in*1e-6
        @variables( urban_city, begin 
        N >= 0.0, (start = L)
        NX >= 0.0, (start = 0.5*L)
        NY >= 0.0, (start = 0.5*L) 
        end)
        η_x = ηx_fixed
        γ_L = γL_fixed
        γ_N = γN_fixed
        ρ_L = ρL_fixed
        ρ_N = ρN_fixed
        R = R_fixed
        I = I_fixed
        T = T_fixed
    end

    # 3. Large city variables
    # Endogenous Variables
    @variables(urban_city, begin
        w, (start = 1.0)        # Wage
        r, (start = 1.0)        # Land Rent
        p, (start = 1.0)        # Price of Y
        x, (start = 1.0); y, (start = 1.0) # Consumption
        X, (start = 1.0); Y, (start = 1.0) # Production Output
        LX, (start = 0.5*L); LY, (start = 0.5*L) # Land Demands
        KX, (start = 0.5); KY, (start = 0.5) # Capital Demands
        K, (start = 1.0)        # Total Capital
    end)
    
    # 4. CES expressions
    # Unit Cost Functions
    @NLexpression(urban_city, 
        cost_X, (γ_L^σ_X * r^(1-σ_X) + γ_N^σ_X * w^(1-σ_X) + (1-γ_L-γ_N)^σ_X * i_bar^(1-σ_X))^(1/(1-σ_X)))
    @NLexpression(urban_city, 
        cost_Y, (ρ_L^σ_Y * r^(1-σ_Y) + ρ_N^σ_Y * w^(1-σ_Y) + (1-ρ_L-ρ_N)^σ_Y * i_bar^(1-σ_Y))^(1/(1-σ_Y)))
    
    # Expenditure Function
    @NLexpression(urban_city, 
        exp_func, (u_bar / Q) * (η_x^σ_D + (1-η_x)^σ_D * p^(1-σ_D))^(1/(1-σ_D)))
    
    # 5. The nonlinear system of equations. These are the MPEC equilbirium constraints.
    # Equilibrium Prices (1b*, 1c*)
    @NLconstraint(urban_city, cost_X / AX == 1.0) 
    @NLconstraint(urban_city, cost_Y == p * AY)
    
    if isnothing(U)
        # Add the Revenue and transfer Constraints
        @NLconstraint(urban_city, R == (r * L) / N)
        @NLconstraint(urban_city, I == (i_bar * K) / N)
        
        # Per-capita tax transfer
        @NLconstraint(urban_city, T == τ * (w + R + I))
    end

    # Consumption Demands (2a*, 2b*)
    @NLconstraint(urban_city, exp_func == (1-τ)*(w + R + I) + T)
    @NLconstraint(urban_city, y == (u_bar/Q) * ((exp_func*Q)/(u_bar*p))^σ_D * (1-η_x)^σ_D)
    @NLconstraint(urban_city, x == (u_bar/Q) * ((exp_func*Q)/(u_bar*1.0))^σ_D * η_x^σ_D)
    
    # Factor Demands (3a*-3f*)
    @NLconstraint(urban_city, NX == (X / AX) * (cost_X / w)^σ_X * γ_N^σ_X)
    @NLconstraint(urban_city, LX == (X / AX) * (cost_X / r)^σ_X * γ_L^σ_X)
    @NLconstraint(urban_city, KX == (X / AX) * (cost_X / i_bar)^σ_X * (1-γ_L-γ_N)^σ_X)
    @NLconstraint(urban_city, NY == (Y / AY) * (cost_Y / w)^σ_Y * ρ_N^σ_Y)
    @NLconstraint(urban_city, LY == (Y / AY) * (cost_Y / r)^σ_Y * ρ_L^σ_Y)
    @NLconstraint(urban_city, KY == (Y / AY) * (cost_Y / i_bar)^σ_Y * (1-ρ_L-ρ_N)^σ_Y)
    
    # Resource Constraints & Market Clearing (4a*-4c*, 6*)
    @NLconstraint(urban_city, N == NX + NY)
    @NLconstraint(urban_city, L == LX + LY)
    @NLconstraint(urban_city, K == KX + KY)
    @NLconstraint(urban_city, Y == N*y)
    
    # 6. Calibrating the large city model
    # Objective: Minimize distance to Table 1 targets
    
    if isnothing(U)
        # Targets from Albouy and Stuart's Table 1
        target_sy = 0.36
        target_θL = 0.025 
        target_θN = 0.825
        target_ϕL = 0.233
        target_ϕN = 0.617
        
        # Uses the square of the sum of the percentage difference
        # (Using just the square differences was failing to match the λs. 
        @NLobjective(urban_city, Min, 
            ( ((p * y) / (x + p * y) - target_sy)/target_sy )^2 +
            ( ((w * NX) / (1.0 * X) - target_θN)/target_θN )^2 +
            ( ((r * LX) / (1.0 * X) - target_θL)/target_θL )^2 +
            ( ((w * NY) / (p * Y)   - target_ϕN)/target_ϕN )^2 +
            ( ((r * LY) / (p * Y)   - target_ϕL)/target_ϕL )^2 +
            ( (NX/N - λ_N) / λ_N)^2 + 
            ( (LX/L - λ_L) / λ_L)^2 
        )
    else
        @NLobjective(urban_city, Min, 0)
    end

    return urban_city
end

model_lg = city_function(N_fixed = 1000.0)
optimize!(model_lg)
termination_status(model_lg) # This tells us where the solver landed
solution_lg = NamedTuple{Tuple(Symbol.(name.(all_variables(model_lg))))}(value.(all_variables(model_lg)))

# 7. Small city comparative statics

support = 0.8:0.01:1.2
results_q = []
model_sm =
        city_function( U = solution_lg.u_bar, Q = 0.8,
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T
        )
set_optimizer_attribute(model_sm, "print_level", 0)

optimize!(model_sm)
push!(results_q, (Q = 0.8, N = value(model_sm[:N]), status = termination_status(model_sm)) )
prev_model = model_sm
for i in support

    global model_sm, prev_model

    model_sm = 
        city_function( U = solution_lg.u_bar, Q = i,
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T
        )
    set_optimizer_attribute(model_sm, "print_level", 0)
    set_start_value.(all_variables(model_sm), value.(all_variables(prev_model)) )
    optimize!(model_sm)

    if termination_status(model_sm) == MOI.LOCALLY_SOLVED || termination_status(model_sm) == MOI.ALMOST_LOCALLY_SOLVED
        prev_model = model_sm
    end
    push!(results_q, (Q = i, N = value(model_sm[:N]), status = termination_status(model_sm)) )

end

plot_x = [row.Q for row in results_q] .- 1.0
plot_y = [row.N for row in results_q] ./ [row.N for row in results_q][plot_x .== 0.0] .- 1.0
Plots.plot(plot_x, plot_y, label = "Q",
            xlims = (-0.2, 0.2), xticks = -0.2:0.05:0.2,
            ylims = (-1.5, 2), yticks = -1.5:0.5:2.0)

results_x = []

model_sm =
        city_function( U = solution_lg.u_bar, AX = 0.8,
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T
        )
set_optimizer_attribute(model_sm, "print_level", 0)

optimize!(model_sm)
push!(results_x, (AX = 0.8, N = value(model_sm[:N]), status = termination_status(model_sm)) )
prev_model = model_sm
for i in support

    global model_sm, prev_model

    model_sm = 
        city_function( U = solution_lg.u_bar, AX = i,
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T
        )
    set_optimizer_attribute(model_sm, "print_level", 0)
    set_start_value.(all_variables(model_sm), value.(all_variables(prev_model)) )
    optimize!(model_sm)

    if termination_status(model_sm) == MOI.LOCALLY_SOLVED || termination_status(model_sm) == MOI.ALMOST_LOCALLY_SOLVED
        prev_model = model_sm
    end
    push!(results_x, (AX = i, N = value(model_sm[:N]), status = termination_status(model_sm)) )

end

plot_x = [row.AX for row in results_x] .- 1.0
plot_y = ( [row.N for row in results_x] ./ [row.N for row in results_x][plot_x .== 0.0] .- 1.0 ) ./ (1-0.36)
Plots.plot(plot_x, plot_y, label = "AX",
            xlims = (-0.2, 0.2), xticks = -0.2:0.05:0.2,
            ylims = (-1.5, 1.5), yticks = -1.5:0.5:1.5
            )

results_y = []

model_sm =
        city_function( U = solution_lg.u_bar, AY = 0.8,
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T
        )
set_optimizer_attribute(model_sm, "print_level", 0)

optimize!(model_sm)
push!(results_y, (AY = 0.8, N = value(model_sm[:N]), status = termination_status(model_sm)) )
prev_model = model_sm
for i in support

    global model_sm, prev_model

    model_sm = 
        city_function( U = solution_lg.u_bar, AY = i,
        ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
        R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T
        )
    set_optimizer_attribute(model_sm, "print_level", 0)
    set_start_value.(all_variables(model_sm), value.(all_variables(prev_model)) )
    optimize!(model_sm)
    
    if termination_status(model_sm) == MOI.LOCALLY_SOLVED || termination_status(model_sm) == MOI.ALMOST_LOCALLY_SOLVED
        prev_model = model_sm
    end
    push!(results_y, (AY = i, N = value(model_sm[:N]), status = termination_status(model_sm)) )

end

plot_x = [row.AY for row in results_y] .- 1.0
plot_y = ( [row.N for row in results_y] ./ [row.N for row in results_y][plot_x .== 0.0] .- 1.0 ) ./ 0.36
Plots.plot(plot_x, plot_y, label = "AY", 
            xlims = (-0.2, 0.2), xticks = -0.2:0.05:0.2,
            ylims = (-2.0, 2.0), yticks = -2.0:0.5:2.0 )
