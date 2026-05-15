using Symbolics

function lower_tri_indices(n)
    indices = Tuple{Int,Int}[]
    for col in 1:n
        for row in col:n
            push!(indices, (row, col))
        end
    end
    return indices
end

function compile_step()

    # model constants/calibration targets. change as needed
    σ_D = 0.667
    σ_X = 0.667
    σ_Y = 0.667
    ι_bar = 1.0
    
    Symbolics.@variables r_s w_s ι_s p_s u_s γL_s γN_s ρL_s ρN_s ηx_s
    
    # cost functions; these can be swapped
    costX_sym = ( γL_s^σ_X * r_s^(1-σ_X) + γN_s^σ_X * w_s^(1-σ_X) + (1 - γL_s - γN_s)^σ_X * ι_bar^( 1-σ_X) )^( 1/(1-σ_X) )
    costY_sym = ( ρL_s^σ_Y * r_s^(1-σ_Y) + ρN_s^σ_Y * w_s^(1-σ_Y) + (1 - ρL_s - ρN_s)^σ_Y * ι_bar^( 1-σ_Y) )^( 1/(1-σ_Y) )
    exp_sym = (ηx_s^σ_D *1.0^(1-σ_D) + (1-ηx_s)^σ_D * p_s^(1-σ_D))^(1/(1-σ_D))

    costX_full = ( γL_s^σ_X * r_s^(1-σ_X) + γN_s^σ_X * w_s^(1-σ_X) + (1 - γL_s - γN_s)^σ_X * ι_s^( 1-σ_X) )^( 1/(1-σ_X) )
    costY_full = ( ρL_s^σ_Y * r_s^(1-σ_Y) + ρN_s^σ_Y * w_s^(1-σ_Y) + (1 - ρL_s - ρN_s)^σ_Y * ι_s^( 1-σ_Y) )^( 1/(1-σ_Y) )
    
    cX_args = [r_s, w_s, γL_s, γN_s]
    cY_args = [r_s, w_s, ρL_s, ρN_s]
    exp_args =[p_s, ηx_s]

    # traded good x, building the function; gradient; jacobian; and hessian

    dcX_sym = Symbolics.gradient(costX_sym, cX_args)
    d2cX_sym = Symbolics.jacobian(dcX_sym, cX_args)

    cX_fun = build_function(costX_sym, cX_args; expression=Val{false})
    dcX_fun_ip! = build_function(dcX_sym, cX_args; expression=Val{false})[2]

    dcX_funs = [build_function(dcX_sym[i], cX_args; expression=Val{false})
                for i in eachindex(dcX_sym)]
    d2cX_row_ip! = [build_function(d2cX_sym[i,:], cX_args; expression=Val{false})[2]
                    for i in axes(d2cX_sym, 1)]
    
    dcX_di_sym = Symbolics.substitute(
                    Symbolics.derivative( costX_full, ι_s ), ι_s => ι_bar )
    d2cX_di_sym = [Symbolics.substitute(
        Symbolics.derivative( Symbolics.derivative(costX_full, ι_s), v),
        ι_s => ι_bar ) for v in cX_args]
    
    dcX_di_fun = build_function( dcX_di_sym, cX_args; expression=Val{false})
    d2cX_di_ip! = build_function( d2cX_di_sym, cX_args; expression=Val{false})[2]
    
    cX_hess_sym = [d2cX_sym[i,j] for (i,j) in lower_tri_indices(4)]
    d2cX_hess_fns = [build_function(cX_hess_sym[i], cX_args; expression = Val{false}) 
                    for i in eachindex(cX_hess_sym)]

    d3cX_r_sym = Symbolics.jacobian( d2cX_sym[1,:], cX_args)
    d3cX_r_hess_sym = [d3cX_r_sym[i,j] for (i,j) in lower_tri_indices(4)]
    d3cX_r_hess_fns = [build_function(d3cX_r_hess_sym[i], cX_args; expression=Val{false}) 
                        for i in eachindex(d3cX_r_hess_sym)]

    d3cX_w_sym = Symbolics.jacobian( d2cX_sym[2,:], cX_args)
    d3cX_w_hess_sym = [d3cX_w_sym[i,j] for (i,j) in lower_tri_indices(4)]
    d3cX_w_hess_fns = [build_function(d3cX_w_hess_sym[i], cX_args; expression=Val{false}) 
                        for i in eachindex(d3cX_w_hess_sym)]
    d3cX_di_sym = Symbolics.jacobian( reshape(d2cX_di_sym, 1, :), cX_args)
    d3cX_di_hess_sym = [d3cX_di_sym[i,j] for (i,j) in lower_tri_indices(4)]
    d3cX_di_hess_fns = [build_function(d3cX_di_hess_sym[i], cX_args; expression=Val{false})
                        for i in eachindex(d3cX_di_hess_sym)]

    # home good y, building the function; gradient; and jacobian
    
    dcY_sym = Symbolics.gradient(costY_sym, cY_args)
    d2cY_sym = Symbolics.jacobian(dcY_sym, cY_args)

    cY_fun = build_function(costY_sym, cY_args; expression=Val{false})
    dcY_fun_ip! = build_function(dcY_sym, cY_args; expression=Val{false})[2]

    dcY_funs = [build_function(dcY_sym[i], cY_args; expression=Val{false})
                for i in eachindex(dcY_sym)]
    d2cY_row_ip! = [build_function(d2cY_sym[i,:], cY_args; expression=Val{false})[2]
                    for i in axes(d2cY_sym, 1)]
    
    dcY_di_sym = Symbolics.substitute(
                    Symbolics.derivative(costY_full, ι_s ), ι_s => ι_bar )
    d2cY_di_sym = [Symbolics.substitute(
        Symbolics.derivative( Symbolics.derivative(costY_full, ι_s), v),
        ι_s => ι_bar ) for v in cY_args]
    
    dcY_di_fun = build_function( dcY_di_sym, cY_args; expression=Val{false})
    d2cY_di_ip! = build_function( d2cY_di_sym, cY_args; expression=Val{false})[2]

    cY_hess_sym = [d2cY_sym[i,j] for (i,j) in lower_tri_indices(4)]
    d2cY_hess_fns = [build_function(cY_hess_sym[i], cY_args; expression = Val{false}) 
                    for i in eachindex(cY_hess_sym)]

    d3cY_r_sym = Symbolics.jacobian( d2cY_sym[1,:], cY_args)
    d3cY_r_hess_sym = [d3cY_r_sym[i,j] for (i,j) in lower_tri_indices(4)]
    d3cY_r_hess_fns = [build_function(d3cY_r_hess_sym[i], cY_args; expression=Val{false}) 
                        for i in eachindex(d3cY_r_hess_sym)]

    d3cY_w_sym = Symbolics.jacobian( d2cY_sym[2,:], cY_args)
    d3cY_w_hess_sym = [d3cY_w_sym[i,j] for (i,j) in lower_tri_indices(4)]
    d3cY_w_hess_fns = [build_function(d3cY_w_hess_sym[i], cY_args; expression=Val{false}) 
                        for i in eachindex(d3cY_w_hess_sym)]

    d3cY_di_sym = Symbolics.jacobian( reshape(d2cY_di_sym, 1, :), cY_args)
    d3cY_di_hess_sym = [d3cY_di_sym[i,j] for (i,j) in lower_tri_indices(4)]
    d3cY_di_hess_fns = [build_function(d3cY_di_hess_sym[i], cY_args; expression=Val{false})
                        for i in eachindex(d3cY_di_hess_sym)]

    # expenditure function, building the function; gradient; and jacobian

    dexp_sym = Symbolics.gradient(exp_sym, exp_args)
    d2exp_sym = Symbolics.jacobian(dexp_sym, exp_args)
    
    exp_fun = build_function(exp_sym, exp_args; expression = Val{false})
    dexp_fun_ip! = build_function(dexp_sym, exp_args; expression = Val{false})[2]

    dexp_dp_fn = build_function( dexp_sym[1], exp_args; expression=Val{false} )
    d2exp_dp_ip! = build_function( d2exp_sym[1,:], exp_args; expression=Val{false} )[2]

    e_hess_sym = [d2exp_sym[i,j] for (i,j) in lower_tri_indices(2)]
    d2e_hess_fns = [build_function(e_hess_sym[i], exp_args; expression=Val{false})
                    for i in eachindex(e_hess_sym)]

    d3e_dp_sym = Symbolics.jacobian( d2exp_sym[1,:], exp_args )
    d3e_dp_hess_sym = [d3e_dp_sym[i,j] for (i,j) in lower_tri_indices(2)]
    d3e_dp_hess_fns = [build_function(d3e_dp_hess_sym[i], exp_args; expression=Val{false})
                        for i in eachindex(d3e_dp_hess_sym)]
    
    return(
        # traded good (X) functions
        cX_fun = cX_fun,
        dcX_fun_ip! = dcX_fun_ip!,
        dcX_funs = dcX_funs, # [1]=∂r, [2]=∂w
        d2cX_row_ip! = d2cX_row_ip!, # [i] fills gradient of ∂cX/∂(arg_i)
        dcX_di_fun = dcX_di_fun,
        d2cX_di_ip! = d2cX_di_ip!,
        # home good (Y) functions
        cY_fun = cY_fun,
        dcY_fun_ip! = dcY_fun_ip!,
        dcY_funs = dcY_funs,
        d2cY_row_ip! = d2cY_row_ip!,
        dcY_di_fun = dcY_di_fun,
        d2cY_di_ip! = d2cY_di_ip!,
        # expenditure functions
        exp_fun = exp_fun,
        dexp_fun_ip! = dexp_fun_ip!,
        dexp_dp_fn = dexp_dp_fn,
        d2exp_dp_ip! = d2exp_dp_ip!,
        # Hessians
        d2cX_hess_fns = d2cX_hess_fns,
        d3cX_r_hess_fns = d3cX_r_hess_fns,
        d3cX_w_hess_fns = d3cX_w_hess_fns,
        d3cX_di_hess_fns = d3cX_di_hess_fns,
        d2cY_hess_fns = d2cY_hess_fns,
        d3cY_r_hess_fns = d3cY_r_hess_fns,
        d3cY_w_hess_fns = d3cY_w_hess_fns,
        d3cY_di_hess_fns = d3cY_di_hess_fns,
        d2e_hess_fns = d2e_hess_fns,
        d3e_dp_hess_fns = d3e_dp_hess_fns
    )

end