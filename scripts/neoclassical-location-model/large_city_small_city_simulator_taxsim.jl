# This version of the script does not include fixed parameters 
# This is the script to call via `include()` when you want to supply your own τ values.
# Assumes you will be loading these packages outside the script; un-comment to change that.
#using JuMP
#using Ipopt
#using Plots

# Using `Ipopt` solver, but others (e.g., KNITRO) could be used if you're licensed.
large_city = Model(Ipopt.Optimizer)
#set_optimizer_attribute(large_city, "print_level", 0) # Un-comment to silence output

# 3. Large city variables
# Endogenous Variables
@variables(large_city, begin
    u_bar, (start = 1.0)    # Utility (Variable in Large City, Fixed in Small)
    w, (start = 1.0)        # Wage
    r, (start = 1.0)        # Land Rent
    p, (start = 1.0)        # Price of Y
    x, (start = 1.0); y, (start = 1.0) # Consumption
    X, (start = 1.0); Y, (start = 1.0) # Production Output
    NX, (start = 500.0); NY, (start = 500.0) # Labor Demands
    LX, (start = 500.0); LY, (start = 500.0) # Land Demands
    KX, (start = 0.5); KY, (start = 0.5) # Capital Demands
    K, (start = 1.0)        # Total Capital
    R, (start = 0.1)
    I, (start = 0.1)
    T, (start = 0.1)
end)

if share_match == "τ"
    @variables(large_city, begin
    0 <= η_x <= 1, (start = 0.5)
    0 <= γ_L <= 1, (start = 0.3)
    0 <= γ_N <= 1, (start = 0.3)
    0 <= ρ_L <= 1, (start = 0.3)
    0 <= ρ_N <= 1, (start = 0.3)
    0 <= τ <= 1, (start = 0.1)
    end)
else
   @variables(large_city, begin
    0 <= η_x <= 1, (start = 0.5)
    0 <= γ_L <= 1, (start = 0.3)
    0 <= γ_N <= 1, (start = 0.3)
    0 <= ρ_L <= 1, (start = 0.3)
    0 <= ρ_N <= 1, (start = 0.3)
    end)
    τ = τ_target
end

# Calibration Parameters (Variables for Large City)


# Exogenous Parameters (Fixed for Large City)
AX = 1.0
AY = 1.0
Q = 1.0

# 4. CES expressions
# Unit Cost Functions
@NLexpression(large_city, 
    cost_X, (γ_L^σ_X * r^(1-σ_X) + γ_N^σ_X * w^(1-σ_X) + (1-γ_L-γ_N)^σ_X * i_bar^(1-σ_X))^(1/(1-σ_X)))
@NLexpression(large_city, 
    cost_Y, (ρ_L^σ_Y * r^(1-σ_Y) + ρ_N^σ_Y * w^(1-σ_Y) + (1-ρ_L-ρ_N)^σ_Y * i_bar^(1-σ_Y))^(1/(1-σ_Y)))

# Expenditure Function
@NLexpression(large_city, 
    exp_func, (u_bar / Q) * (η_x^σ_D + (1-η_x)^σ_D * p^(1-σ_D))^(1/(1-σ_D)))

# 5. The nonlinear system of equations. These are the MPEC equilbirium constraints.
# Equilibrium Prices (1b*, 1c*)
@NLconstraint(large_city, cost_X / AX == 1.0) 
@NLconstraint(large_city, cost_Y == p * AY)


# Add the Revenue and transfer Constraints
@NLconstraint(large_city, R == (r * L_large) / N_large)
@NLconstraint(large_city, I == (i_bar * K) / N_large)

# Per-capita tax transfer
@NLconstraint(large_city, T == τ * (w + R + I))

# Consumption Demands (2a*, 2b*)
@NLconstraint(large_city, exp_func == (1-τ)*(w + R + I) + T)
@NLconstraint(large_city, y == (u_bar/Q) * ((exp_func*Q)/(u_bar*p))^σ_D * (1-η_x)^σ_D)
@NLconstraint(large_city, x == (u_bar/Q) * ((exp_func*Q)/(u_bar*1.0))^σ_D * η_x^σ_D)

# Factor Demands (3a*-3f*)
@NLconstraint(large_city, NX == (X / AX) * (cost_X / w)^σ_X * γ_N^σ_X)
@NLconstraint(large_city, LX == (X / AX) * (cost_X / r)^σ_X * γ_L^σ_X)
@NLconstraint(large_city, KX == (X / AX) * (cost_X / i_bar)^σ_X * (1-γ_L-γ_N)^σ_X)
@NLconstraint(large_city, NY == (Y / AY) * (cost_Y / w)^σ_Y * ρ_N^σ_Y)
@NLconstraint(large_city, LY == (Y / AY) * (cost_Y / r)^σ_Y * ρ_L^σ_Y)
@NLconstraint(large_city, KY == (Y / AY) * (cost_Y / i_bar)^σ_Y * (1-ρ_L-ρ_N)^σ_Y)

# Resource Constraints & Market Clearing (4a*-4c*, 6*)
@NLconstraint(large_city, N_large == NX + NY)
@NLconstraint(large_city, L_large == LX + LY)
@NLconstraint(large_city, K == KX + KY)
@NLconstraint(large_city, Y == N_large*y)

# 6. Calibrating the large city model
# Objective: Minimize distance to Table 1 targets

# Uses the square of the sum of the percentage difference
# (Using just the square differences was failing to match the λs.) 

if share_match == "τ" # This breaks if you try to supply τ along with λs... Which makes sense!
    @NLobjective(large_city, Min, 
        ( ((p * y) / (x + p * y) - target_sy)/target_sy )^2 +
        ( ((w * NX) / (1.0 * X) - target_θN)/target_θN )^2 +
        ( ((r * LX) / (1.0 * X) - target_θL)/target_θL )^2 +
        ( ((w * NY) / (p * Y)   - target_ϕN)/target_ϕN )^2 +
        ( ((r * LY) / (p * Y)   - target_ϕL)/target_ϕL )^2 +
        ( (τ - τ_target) / τ_target)^2
    )
else
    @NLobjective(large_city, Min, 
        ( ((p * y) / (x + p * y) - target_sy)/target_sy )^2 +
        ( ((w * NX) / (1.0 * X) - target_θN)/target_θN )^2 +
        ( ((r * LX) / (1.0 * X) - target_θL)/target_θL )^2 +
        ( ((w * NY) / (p * Y)   - target_ϕN)/target_ϕN )^2 +
        ( ((r * LY) / (p * Y)   - target_ϕL)/target_ϕL )^2 +
        ( (NX/N_large - λ_N) / λ_N)^2 + 
        ( (LX/L_large - λ_L) / λ_L)^2 
    )
end
optimize!(large_city)
status = termination_status(large_city) # This tells us where the solver landed

if status == MOI.LOCALLY_SOLVED || status == MOI.ALMOST_LOCALLY_SOLVED
    # Store these as standard global variables first
    η_val = value(η_x)
    γL_val = value(γ_L)
    γN_val = value(γ_N)
    ρL_val = value(ρ_L)
    ρN_val = value(ρ_N)
    u_val = value(u_bar)
    R_val = value(R)
    I_val = value(I)
    T_val = value(T)
    if share_match == "τ"
        τ_val = value(τ)
        println("Calibration Successful. τ_val = ", τ_val)
    else
        τ_val = τ_target
        println("Calibration Successful. η_x = ", η_val)
    end
else
    error("Calibration failed with status: ", status)
end
# 7. Small city comparative statics

# Initialize storage bin for all the completed paramter "experiments"
meta_results = []

for var in ["Q", "AX", "AY"]
    # Initialize storage bin to collect values
    results = []
    
    for val in 0.8:0.01:1.2

        small_city = Model(Ipopt.Optimizer)
        set_optimizer_attribute(small_city, "print_level", 0)

        # We have to put _small after a lot of these, or Julia throws us a lot of scope warnings.
        # They don't affect anything in this script, but are annoying, and could concievably be an issue if used in a bigger project.
        if var == "Q"
            AX_small = 1.0
            AY_small = 1.0
            Q_small = val
        elseif var == "AX"
            AX_small = val
            AY_small = 1.0
            Q_small = 1.0
        elseif var == "AY"
            AX_small = 1.0
            AY_small = val
            Q_small = 1.0
        else
            error("Please supply 'AX,' 'AY,' or 'Q' as arguments")
        end
    
        # Initial values helpfully provided by Stuart! They seem to work very well here, too.
        if var == "Q" && 0.8 <= val <= 0.84
            @variables(small_city, begin
                N_small, (start = 0.0001) # Total Population
                w_small, (start = 0.65)        # Wage
                p_small, (start = 0.65)        # Price of Y
                r_small, (start = 0.005)       # Land rent
                x_small, (start = 0.55); y_small, (start = 0.55) # Consumption
                X_small, (start = 0.00001); Y_small, (start = 0.00001) # Production Output
                NX_small, (start = 0.00001); NY_small, (start = 0.00001) # Labor Demands
                LX_small, (start = 0.0004); LY_small, (start = 0.0004) # Land Demands
                KX_small, (start = 0.000001); KY_small, (start = 0.000001) # Capital Demands
                K_small, (start = 0.000001)        # Total Capital
            end)
        elseif var == "Q" && 0.841 <= val <= 1.0
            @variables(small_city, begin
                N_small, (start = 0.0001) # Total Population
                w_small, (start = 0.7)        # Wage
                p_small, (start = 0.7)        # Price of Y
                r_small, (start = 0.05)       # Land rent
                x_small, (start = 0.4); y_small, (start = 0.4) # Consumption
                X_small, (start = 0.0003); Y_small, (start = 0.0003) # Production Output
                NX_small, (start = 0.0003); NY_small, (start = 0.0003) # Labor Demands
                LX_small, (start = 0.0004); LY_small, (start = 0.0004) # Land Demands
                KX_small, (start = 0.00001); KY_small, (start = 0.00001) # Capital Demands
                K_small, (start = 0.00001)        # Total Capital
            end)
        
        elseif var == "Q" && 1.001 <= val <= 1.2
            @variables(small_city, begin
                N_small, (start = 0.001) # Total Population
                w_small, (start = 0.8)        # Wage
                p_small, (start = 0.8)        # Price of Y
                r_small, (start = 0.07)       # Land rent
                x_small, (start = 0.4); y_small, (start = 0.4) # Consumption
                X_small, (start = 0.0001); Y_small, (start = 0.0001) # Production Output
                NX_small, (start = 0.0001); NY_small, (start = 0.0001) # Labor Demands
                LX_small, (start = 0.0001); LY_small, (start = 0.0001) # Land Demands
                KX_small, (start = 0.0001); KY_small, (start = 0.0001) # Capital Demands
                K_small, (start = 0.0001)        # Total Capital
            end)
        
        elseif var == "AX" && 0.8 <= val <= 0.95
            @variables(small_city, begin
                N_small, (start = 0.0001) # Total Population
                w_small, (start = 0.5)        # Wage
                p_small, (start = 0.5)        # Price of Y
                r_small, (start = 0.01)       # Land rent
                x_small, (start = 0.5); y_small, (start = 0.5) # Consumption
                X_small, (start = 0.0001); Y_small, (start = 0.0001) # Production Output
                NX_small, (start = 0.0001); NY_small, (start = 0.0001) # Labor Demands
                LX_small, (start = 0.0004); LY_small, (start = 0.0004) # Land Demands
                KX_small, (start = 0.00001); KY_small, (start = 0.00001) # Capital Demands
                K_small, (start = 0.00001)        # Total Capital
            end)
        
        elseif var == "AX" && 0.951 <= val <= 1.15
            @variables(small_city, begin
                N_small, (start = 0.0001) # Total Population
                w_small, (start = 0.7)        # Wage
                p_small, (start = 0.7)        # Price of Y
                r_small, (start = 0.05)       # Land rent
                x_small, (start = 0.5); y_small, (start = 0.5) # Consumption
                X_small, (start = 0.0001); Y_small, (start = 0.0001) # Production Output
                NX_small, (start = 0.0001); NY_small, (start = 0.0001) # Labor Demands
                LX_small, (start = 0.0004); LY_small, (start = 0.0004) # Land Demands
                KX_small, (start = 0.0001); KY_small, (start = 0.0001) # Capital Demands
                K_small, (start = 0.0001)        # Total Capital
            end)
            
        elseif var == "AX" && 1.151 <= val <= 1.19
            @variables(small_city, begin
                N_small, (start = 0.001) # Total Population
                w_small, (start = 0.8)        # Wage
                p_small, (start = 1.2)        # Price of Y
                r_small, (start = 0.1)       # Land rent
                x_small, (start = 0.45); y_small, (start = 0.45) # Consumption
                X_small, (start = 0.0008); Y_small, (start = 0.0008) # Production Output
                NX_small, (start = 0.0008); NY_small, (start = 0.0008) # Labor Demands
                LX_small, (start = 0.0005); LY_small, (start = 0.0005) # Land Demands
                KX_small, (start = 0.0001); KY_small, (start = 0.0001) # Capital Demands
                K_small, (start = 0.0001)        # Total Capital
            end)
        
        elseif var == "AX" && val == 1.2
            @variables(small_city, begin
                N_small, (start = 0.001) # Total Population
                w_small, (start = 0.8)        # Wage
                p_small, (start = 1.3)        # Price of Y
                r_small, (start = 0.2)       # Land rent
                x_small, (start = 0.6); y_small, (start = 0.3) # Consumption
                X_small, (start = 0.001); Y_small, (start = 0.001) # Production Output
                NX_small, (start = 0.0008); NY_small, (start = 0.0008) # Labor Demands
                LX_small, (start = 0.0005); LY_small, (start = 0.0005) # Land Demands
                KX_small, (start = 0.0001); KY_small, (start = 0.0001) # Capital Demands
                K_small, (start = 0.0001)        # Total Capital
            end)
        
        elseif var == "AY" && 0.8 <= val <= 0.92
            @variables(small_city, begin
                N_small, (start = 0.0001) # Total Population
                w_small, (start = 0.7)        # Wage
                p_small, (start = 0.7)        # Price of Y
                r_small, (start = 0.005)       # Land rent
                x_small, (start = 0.55); y_small, (start = 0.55) # Consumption
                X_small, (start = 0.0001); Y_small, (start = 0.0001) # Production Output
                NX_small, (start = 0.0001); NY_small, (start = 0.0001) # Labor Demands
                LX_small, (start = 0.0004); LY_small, (start = 0.0004) # Land Demands
                KX_small, (start = 0.00001); KY_small, (start = 0.00001) # Capital Demands
                K_small, (start = 0.00001)        # Total Capital
            end)
        
        elseif var == "AY" && 0.921 <= val <= 1.03
            @variables(small_city, begin
                N_small, (start = 0.0006) # Total Population
                w_small, (start = 0.6)        # Wage
                p_small, (start = 0.8)        # Price of Y
                r_small, (start = 0.05)       # Land rent
                x_small, (start = 0.45); y_small, (start = 0.45) # Consumption
                X_small, (start = 0.0003); Y_small, (start = 0.0003) # Production Output
                NX_small, (start = 0.0003); NY_small, (start = 0.0003) # Labor Demands
                LX_small, (start = 0.0004); LY_small, (start = 0.0004) # Land Demands
                KX_small, (start = 0.0001); KY_small, (start = 0.00001) # Capital Demands
                K_small, (start = 0.0001)        # Total Capital
            end)
        else 
            @variables(small_city, begin
                N_small, (start = 0.0001) # Total Population
                w_small, (start = 0.6)        # Wage
                p_small, (start = 0.8)        # Price of Y
                r_small, (start = 0.05)       # Land rent
                x_small, (start = 0.45); y_small, (start = 0.45) # Consumption
                X_small, (start = 0.0003); Y_small, (start = 0.0003) # Production Output
                NX_small, (start = 0.0003); NY_small, (start = 0.0003) # Labor Demands
                LX_small, (start = 0.0004); LY_small, (start = 0.0004) # Land Demands
                KX_small, (start = 0.0001); KY_small, (start = 0.0001) # Capital Demands
                K_small, (start = 0.0001)        # Total Capital
            end)
        end

        # CES Expressions for small city
        # Unit Cost Functions
        @NLexpression(small_city, cost_X_small, (γL_val^σ_X * r_small^(1-σ_X) + γN_val^σ_X * w_small^(1-σ_X) + (1-γL_val-γN_val)^σ_X * i_bar^(1-σ_X))^(1/(1-σ_X)))
        @NLexpression(small_city, cost_Y_small, (ρL_val^σ_Y * r_small^(1-σ_Y) + ρN_val^σ_Y * w_small^(1-σ_Y) + (1-ρL_val-ρN_val)^σ_Y * i_bar^(1-σ_Y))^(1/(1-σ_Y)))
    
        # Expenditure Function
        @NLexpression(small_city, exp_func_small, (u_val / Q_small) * (η_val^σ_D + ( (1-η_val)^σ_D * p_small^(1-σ_D) ) )^(1/(1-σ_D)))

        # Small city constraints
        # Equilibrium Prices
        @NLconstraint(small_city, cost_X_small / AX_small == 1.0) 
        @NLconstraint(small_city, cost_Y_small == p_small * AY_small)
    
        # Consumption Demands
        @NLconstraint(small_city, exp_func_small == (1-τ_val)*(w_small + R_val + I_val) + T_val)
        @NLconstraint(small_city, y_small == (u_val/Q_small) * ((exp_func_small*Q_small)/(u_val*p_small))^σ_D * (1-η_val)^σ_D)
        @NLconstraint(small_city, x_small == (u_val/Q_small) * ((exp_func_small*Q_small)/(u_val*1.0))^σ_D * η_val^σ_D)
    
        # Factor Demands
        @NLconstraint(small_city, NX_small == ( X_small / AX_small ) * ( cost_X_small / w_small )^σ_X * γN_val^σ_X)
        @NLconstraint(small_city, LX_small == ( X_small / AX_small ) * ( cost_X_small / r_small )^σ_X * γL_val^σ_X)
        @NLconstraint(small_city, KX_small == ( X_small / AX_small ) * ( cost_X_small / i_bar )^σ_X * (1-γL_val-γN_val)^σ_X)
        @NLconstraint(small_city, NY_small == ( Y_small / AY_small ) * ( cost_Y_small / w_small )^σ_Y * ρN_val^σ_Y)
        @NLconstraint(small_city, LY_small == ( Y_small / AY_small ) * ( cost_Y_small / r_small )^σ_Y * ρL_val^σ_Y)
        @NLconstraint(small_city, KY_small == ( Y_small / AY_small ) * ( cost_Y_small / i_bar )^σ_Y * (1-ρL_val-ρN_val)^σ_Y)
    
        # Resource Constraints & Market Clearing
        @NLconstraint(small_city, N_small == NX_small + NY_small)
        @NLconstraint(small_city, L_small == LX_small + LY_small) # L_large fixed at 1000 for Large City
        @NLconstraint(small_city, K_small == KX_small + KY_small)
        @NLconstraint(small_city, Y_small == N_small*y_small)

        # No need to optimize over an objective because we are supplying all the parameters now.
        optimize!(small_city)
        stat_small = termination_status(small_city)
    
        if var == "Q"
            push!( results, (Q=val, N=value(N_small), solved=stat_small ) )
        elseif var == "AX"
            push!( results, (AX=val, N=value(N_small), solved=stat_small ) )
        else
            push!( results, (AY=val, N=value(N_small), solved=stat_small ) )
        end

    end

    # Extract values into a vector

    if var == "Q"
        output_matrix = 
            [ [row.Q for row in results] [row.N for row in results] [row.solved for row in results] ]
    elseif var == "AX"
        output_matrix = 
            [ [row.AX for row in results] [row.N for row in results] [row.solved for row in results] ]
        x_lab = "Percent change in AX"
    else 
        output_matrix = 
            [ [row.AY for row in results] [row.N for row in results] [row.solved for row in results] ]
        x_lab = "Percent change in AY"
    end

    output_matrix[:,2] = output_matrix[:,2]./output_matrix[output_matrix[:,1] .== 1.0, 2] .- 1.0
    output_matrix[:,1] = output_matrix[:,1] .- 1.0

    push!(meta_results, output_matrix)
end

# Produce Figure A.3 from Albouy and Stuart
p1 = 
    Plots.plot(meta_results[1][:,1], meta_results[1][:,2], 
    xlabel = "Percent change in Q", ylabel = "Percent change in N", legend = false, 
    title = "Quality of life:")
p2 = 
    Plots.plot(meta_results[2][:,1], 
    meta_results[2][:,2]./ (1 - target_sy), # Based on footnote 33, I believe this normaization is correct.
    xlabel = "Percent change in AX", ylabel = "Percent change in N", legend = false, 
    title = "Trade productivity:")
p3 = 
    Plots.plot(meta_results[3][:,1], meta_results[3][:,2]./target_sy, # See comment on p2
    xlabel = "Percent change in AY", ylabel = "Percent change in N", legend = false, 
    title = "Home productivity:")

dw, dh = default(:size)
small_plot = plot(p1,p2,p3, layout = (3,1), legend = false, size = (dw, dh * 3))
