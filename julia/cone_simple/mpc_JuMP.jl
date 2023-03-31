function mpc_JuMP(optimizer, params, x0, A, B, f)
    Nh = params.N
    nx = params.nx
    nu = params.nu
    α_max = params.c_cone[3]
    NN = Nh*nx + (Nh-1)*nu
    
    inds = reshape(1:(nx+nu)*Nh,nx+nu,Nh)  
    xinds = [z[1:nx] for z in eachcol(inds)]
    uinds = [z[nx+1:end] for z in eachcol(inds)][1:Nh-1]    
    
    model = Model(optimizer)
    
    @variable(model, z[1:NN])  # z is all decision variables (X U)
    
    P = zeros(NN, NN)
    q = zeros(NN, 1) 
    # Cost function   
    for j = 1:Nh-1
        P[(j-1)*(nx+nu).+(1:nx),(j-1)*(nx+nu).+(1:nx)], q[(j-1)*(nx+nu).+(1:nx)], 
        P[(j)*(nx+nu).+(1:nu),(j)*(nx+nu).+(1:nu)], q[(j)*(nx+nu).+(1:nu)] = stage_cost_expansion(params, j)
    end    
    P[end-nx+1:end,end-nx+1:end], q[end-nx+1:end] = term_cost_expansion(params)
    @objective(model, Min, 0.5*dot(z,P,z) + dot(q,z))
    
    # Dynamics Constraints
    for k = 1:Nh-1
        @constraint(model, A*z[xinds[k]] .+ B*z[uinds[k]] .+ f .== z[xinds[k+1]])
    end
    
    # Initial condition 
    @constraint(model, z[xinds[1]] .== x0)
    
    # Thrust angle constraint
    for k = 1:Nh-1
        u1,u2,u3 = z[uinds[k]]
        @constraint(model, [α_max * u3, u1, u2] in JuMP.SecondOrderCone())
    end
    
    # # Goal constraint
    # if goal_constraint 
    #     @constraint(model, z[xinds[N]] .== prob.xf)
    # end    
    
    optimize!(model)   
    termination_status(model) == INFEASIBLE && print("Other solver says INFEASIBLE\n")
    return value.(z[uinds[1]])[:]
end