using JuMP, Ipopt
using Symbolics

# load up our functions
# note: they're separated so you can tweak the "_func_compiler.jl" with new functions, if desired
include("symbolic_function_city_simulator_func_compiler.jl")
include("symbolic_function_city_simulator_builder.jl")

# loads the compiled cost/expenditure functions; outputting `fns`
fns = compile_step()

# takes `fns` from `compile_step()` and builds the JuMP model
model_lg = build_model(fns;
    Q = 1.0, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
    τ = 0.0, ι_bar = 1.0, 
    N_fixed = 1000.0
)

# solve the model output and acquire the solutions
optimize!(model_lg)

solution_mat_lg = [all_variables(model_lg) value.(all_variables(model_lg))]
solution_lg = NamedTuple{Tuple(Symbol.(name.(all_variables(model_lg))))}(value.(all_variables(model_lg)))
@show termination_status(model_lg)
@show objective_value(model_lg)

# diagnostics
inc = solution_lg.x + solution_lg.p * solution_lg.y
@show solution_lg.γ_L solution_lg.γ_N solution_lg.ρ_L solution_lg.ρ_N solution_lg.η_x
@show solution_lg.p * solution_lg.y / inc           # target: 0.36
@show solution_lg.r * solution_lg.LX / solution_lg.X      # target: 0.025
@show solution_lg.w * solution_lg.NX / solution_lg.X      # target: 0.825
@show solution_lg.r * solution_lg.LY / (solution_lg.p * solution_lg.Y)  # target: 0.233
@show solution_lg.w * solution_lg.NY / (solution_lg.p * solution_lg.Y)  # target: 0.617

# input solution values from large city; use to solve small city values
model_sm = build_model(fns;
    Q = 1.0, A_X = 1.0, A_Y = 1.0, L_0 = 1000.0, 
    τ = 0.0, ι_bar = 1.0, 
    ηx_fixed = solution_lg.η_x, γL_fixed = solution_lg.γ_L, γN_fixed = solution_lg.γ_N, ρL_fixed = solution_lg.ρ_L, ρN_fixed = solution_lg.ρ_N,
    R_fixed = solution_lg.R, I_fixed = solution_lg.I, T_fixed = solution_lg.T,
    u_bar = solution_lg.u_bar_v,
)

optimize!(model_sm)

solution_mat_sm = [all_variables(model_sm) value.(all_variables(model_sm))]
solution_sm = NamedTuple{Tuple(Symbol.(name.(all_variables(model_lg))))}(value.(all_variables(model_lg)))
@show termination_status(model_sm)
@show objective_value(model_sm)
