function stage_cost(p::NamedTuple, x, u, k)
    dx = x - p.Xref[k]
    du = u - p.Uref[k]
    return 0.5 * dx' * p.Q * dx + 0.5 * du' * p.R * du
end
function term_cost(p::NamedTuple, x)
    dx = x - p.Xref[p.N]
    return 0.5 * dx' * p.Qf * dx
end
function stage_cost_expansion(p::NamedTuple, k)
    dx = -p.Xref[k]
    du = -p.Uref[k]
    return p.Q, p.Q * dx, p.R, p.R * du  # Hessian and gradient
end
function term_cost_expansion(p::NamedTuple)
    dx = -p.Xref[p.N]
    return p.Qf, p.Qf * dx
end
function backward_pass!(params, X, U, P, p, d, K, reg, μ, μx, ρ, λ, λc)
    """LQR backward pass with AL. This !function update its parameters
    """
    N = params.N
    ΔJ = 0.0    # expected cost reduction

    # terminal cost expansion
    P[N], p[N] = term_cost_expansion(params)

    if (params.ncx > 0)
        # # add AL terms for the state constraint at the final time step
        hxv = ineq_con_x(params, X[N])  # h(x) violation  
        mask = eval_mask(μx[N], hxv)
        ∇hx = ineq_con_x_jac(params, X[N])
        # add these into the cost-to-go p and P
        p[N] += ∇hx' * (μx[N] - ρ * (mask * [params.x_max; -params.x_min]))  # multiplier term (1st)
        P[N] += ρ * ∇hx' * mask * ∇hx                # penalty term (2nd)
    end

    if (params.ncg > 0)
        # add AL terms for goal constraint 
        ∇hx = diagm(ones(params.nx))
        # add these into the CTG p and P (equality active)
        p[N] += ∇hx' * (λ - ρ * params.Xref[N])
        P[N] += ρ * ∇hx'∇hx
    end

    # iterate from N-1 to 1 backwards
    for k = (N-1):(-1):1
        # dynamics jacobians (linear) (still give the same A, B, f)
        A = ForwardDiff.jacobian(_x -> discrete_dynamics(params, _x, U[k], k), X[k])
        B = ForwardDiff.jacobian(_u -> discrete_dynamics(params, X[k], _u, k), U[k])
        f = discrete_dynamics(params, zeros(params.nx), zeros(params.nu), k)

        Sxx, Sx, Suu, Su = stage_cost_expansion(params, k)

        # one-step cost expansion: Q, Q*dx, R, R*du
        Sx += A' * (P[k+1] * f + p[k+1])
        Su += B' * (P[k+1] * f + p[k+1])
        Sxx += A' * (P[k+1]) * A
        Suu += B' * (P[k+1] + reg * I) * B
        Sux = B' * (P[k+1]) * A

        if (params.ncu > 0)
            # control constraints
            huv = ineq_con_u(params, U[k])  # calculate h(u) constraint
            mask = eval_mask(μ[k], huv)  # choose active
            ∇hu = ineq_con_u_jac(params, U[k])
            Su += ∇hu' * (μ[k] - ρ * (mask * [params.u_max; -params.u_min])) # add to cost
            Suu += ρ * ∇hu' * mask * ∇hu
        end
        if (params.ncx > 0)
            # state constraints
            hxv = ineq_con_x(params, X[k])
            mask = eval_mask(μx[k], hxv)
            ∇hx = ineq_con_x_jac(params, X[k])
            Sx += ∇hx' * (μx[k] - ρ * (mask * [params.x_max; -params.x_min]))
            Sxx += ρ * ∇hx' * mask * ∇hx
        end
        if (params.ncu_cone > 0)
            # conic constraints
            Qu, Quu = conic_cost_expansion(params, U[k], λc[k], ρ * cone_scale, k)
            # display(Quu)
            Su += Qu
            Suu += Quu
        end

        F = cholesky(Symmetric(Suu))
        d[k] = F \ Su
        K[k] = F \ Sux

        # Cost-to-go Recurrence (PSD stabilizing version, last term)
        P[k] = Sxx + K[k]' * Suu * K[k] - K[k]' * Sux - Sux' * K[k]
        p[k] = Sx + K[k]' * Suu * d[k] - K[k]' * Su - Sux' * d[k]
    end

    return 0.0 
end
function trajectory_AL_cost(params, X, U, μ, μx, ρ, λ, λc)
    # Evaluate merit function
    N = params.N
    J = 0.0
    for k = 1:N-1
        J += stage_cost(params, X[k], U[k], k)
        if params.ncu > 0
            # AL terms for ineq_con_u
            huv = ineq_con_u(params, U[k])
            mask = eval_mask(μ[k], huv)
            J += dot(μ[k], huv) + 0.5 * ρ * huv' * mask * huv
        end
        if params.ncx > 0
            # AL terms for ineq_con_x
            hxv = ineq_con_x(params, X[k])
            mask = eval_mask(μx[k], hxv)
            J += dot(μx[k], hxv) + 0.5 * ρ * hxv' * mask * hxv
        end
    end
    if params.ncx > 0
        # AL terms for state constraint at last time step
        J += term_cost(params, X[N])
        hxv = ineq_con_x(params, X[params.N])
        mask = eval_mask(μx[params.N], hxv)
        J += dot(μx[params.N], hxv) + 0.5 * ρ * hxv' * mask * hxv
    end
    if params.ncg > 0
        # AL terms for goal constraint
        hxv = X[N] - params.Xref[N]
        J += dot(λ, hxv) + 0.5 * ρ * hxv' * hxv
    end
    return J
end
function forward_pass!(params, X, U, K, d, ΔJ, Xn, Un, μ, μx, ρ, λ, λc; max_linesearch_iters=20)
    """LQR forward pass
    This !function update its parameters
    """
    N = params.N
    # Forward Rollout
    Xn[1] .= X[1]
    J = trajectory_AL_cost(params, X, U, μ, μx, ρ, λ, λc)
    for k = 1:(N-1)
        Un[k] = -d[k] - K[k] * Xn[k]
        Xn[k+1] = discrete_dynamics(params, Xn[k], Un[k], k)
    end
    Jn = trajectory_AL_cost(params, Xn, Un, μ, μx, ρ, λ, λc)
    ΔJ = J - Jn  # Just experiment, not the right way
    X .= Xn
    U .= Un
    return Jn, ΔJ
end
function eval_mask(μv, huv)
    # Extract active inequality constraints
    # active set mask
    mask = Diagonal(zeros(length(huv)))
    for i = 1:length(huv)
        mask[i, i] = (μv[i] > 0 || huv[i] > 0)
    end
    mask
end
function solve!(params, X, U, P, p, K, d, Xn, Un; atol=1e-3, max_iters=50, max_inner_iters=10, verbose=true, ρ=1, ρ_max=1e8, ϕ=10)

    # first check the sizes of everything
    # @assert length(X) == params.N
    # @assert length(U) == params.N-1
    # @assert length(X[1]) == params.nx
    # @assert length(U[1]) == params.nu
    ΔJ = 0.0
    # initial rollout
    N = params.N
    for i = 1:N-1
        X[i+1] = discrete_dynamics(params, X[i], U[i], i)
    end

    reg_min = 1e-8  # can be zero
    reg = reg_min

    μ = params.μ
    μx = params.μx
    λ = params.λ
    λc = params.λc

    # Solve unconstrained LQR first
    backward_pass!(params, X, U, P, p, d, K, reg, μ, μx, 0.0, λ, λc)
    J, ΔJ = forward_pass!(params, X, U, K, d, ΔJ, Xn, Un, μ, μx, 0.0, λ, λc)

    # Inner loop to solve unconstrained problem (Riccati)
    for iter = 1:max_iters
        α = 1.0  # no line-search

        # Riccati solve
        for i = 1:max_inner_iters
            backward_pass!(params, X, U, P, p, d, K, reg, μ, μx, ρ, λ, λc)
            J, ΔJ = forward_pass!(params, X, U, K, d, ΔJ, Xn, Un, μ, μx, ρ, λ, λc)
            if verbose 
                @show ΔJ
            end
            if 0.0 <= ΔJ <= 1.0  # experiment, not a good one
                break
            end
            # reg = reg*10.0  # if needed
        end
        reg = reg_min

        if verbose
            if rem(iter - 1, 10) == 0
                @printf "iter     J           ΔJ        |d|         α        reg         ρ\n"
                @printf "---------------------------------------------------------------------\n"
            end
            @printf("%3d   %10.3e  %9.2e  %9.2e  %6.4f   %9.2e   %9.2e\n",
                iter, J, ΔJ, 0, α, reg, ρ)
        end

        # AL update and check constraint violation
        convio = 0.0

        if (params.ncu > 0)
            # control constraints (inequality)
            for k = 1:N-1
                huv = ineq_con_u(params, U[k])
                mask = eval_mask(μ[k], huv)
                # update dual
                convio = max(convio, norm(huv + abs.(huv), Inf))
                μ[k] .= max.(0, μ[k] - ρ * (mask * [params.u_max; -params.u_min]))
            end
        end
        if (params.ncx > 0)
            # state constraints (inequality)
            for k = 1:N
                hxv = ineq_con_x(params, X[k])
                mask = eval_mask(μx[k], hxv)
                # update dual
                convio = max(convio, norm(hxv + abs.(hxv), Inf))
                μx[k] .= max.(0, μx[k] - ρ * (mask * [params.x_max; -params.x_min]))
            end
        end
        if (params.ncu_cone > 0)
            # conic constraints 
            if verbose
                print("update cone\n")
            end
            for k = 1:N-1
                Uc = cone_u(params, U[k])
                huc = norm(Uc[1:2]) - Uc[3]
                convio = max(convio, norm(huc + abs.(huc), Inf))
                # update dual
                λc[k] .= projection(λc[k] - Uc * cone_scale * ρ)
            end
        end
        if (params.ncg > 0)
            # goal constraint (equality)
            hxv = X[N] - params.Xref[N]
            λ .-= ρ * params.Xref[N]
            convio = max(convio, norm(hxv, Inf))
        end

        ρ *= ϕ  # update penalty

        if verbose
            @show convio
        end
        if convio < atol  # if terminal condition with contraint violation
            # print("Our solver says SUCCESS\n") # @info "success!"
            return Un[1]
        end
        if ρ > ρ_max
            print("MAX PENALTY\n")
            return Un[1]
        end
    end
    print("MAX ITER\n")
    return Un[1]  # get full U and X by from passed argument
end

# =============================
# Constraints
# =============================
function ineq_con_x(p, x)
    [x - p.x_max; -x + p.x_min]
end
function ineq_con_u(p, u)
    [u - p.u_max; -u + p.u_min]
end
function ineq_con_u_jac(params, u)
    ForwardDiff.jacobian(_u -> ineq_con_u(params, _u), u)  # lazy way
end
function ineq_con_x_jac(p, x)
    ForwardDiff.jacobian(_x -> ineq_con_x(p, _x), x)  # lazy way
end
function mat_from_vec(X::Vector{Vector{Float64}})::Matrix
    # convert a vector of vectors to a matrix 
    Xm = hcat(X...)
    return Xm
end
function shift_fill(U::Vector)
    N = length(U)
    for k = 1:N-1
        U[k] .= U[k+1]
    end
    U[N] .= U[N-1]
end
function discrete_dynamics(params::NamedTuple,x,u,k)
    A = [1.0 0.0 0.0 0.05 0.0 0.0; 
        0.0 1.0 0.0 0.0 0.05 0.0; 
        0.0 0.0 1.0 0.0 0.0 0.05; 
        0.0 0.0 0.0 1.0 0.0 0.0; 
        0.0 0.0 0.0 0.0 1.0 0.0; 
        0.0 0.0 0.0 0.0 0.0 1.0]
    B = [0.000125 0.0 0.0; 
        0.0 0.000125 0.0;
        0.0 0.0 0.000125; 
        0.005 0.0 0.0; 
        0.0 0.005 0.0; 
        0.0 0.0 0.005]
    f = [0.0, 0.0, -0.0122625, 0.0, 0.0, -0.4905]
    return A*x + B*u + f
end

function main()
    nx = 6
    nu = 3
    N = 101
    dt = 0.1
    t_vec = dt*(0:N-1)
    x0 = [4, 2, 20, -3, 2, -5.0]
    xg = [0,0,0,0,0,0.0]
    Xref = [deepcopy(xg) for i = 1:N]
    Uref = [zeros(nu) for i = 1:N-1]

    Q = 1e-2*Diagonal([1,1,1,1.0,1,1])
    R = 1e-1*Diagonal([1,1,1])
    Qf = 1000*Q

    u_min = -170*ones(nu)
    u_max =  170*ones(nu)

    # state is x y v θ
    x_min = [-2,-2,-2,-2]
    x_max = [6,8,3,2]

    ncx = 2*nx*0
    ncu = 2*nu
    ncg = 1
    ncu_cone = nu; 

    params = (
        nx = nx,
        nu = nu,
        ncx = ncx,
        ncu = ncu,
        ncg = ncg,
        ncu_cone = ncu_cone,
        A_cone = A_cone,
        c_cone = c_cone,
        N = N,
        Q = Q,
        R = R,
        Qf = Qf,
        u_min = u_min,
        u_max = u_max,
        x_min = x_min,
        x_max = x_max,
        Xref = Xref,
        Uref = Uref,
        dt = dt,
        mc = 1.0,
        mp = 0.2,
        l = 0.5,
        g = 9.81,
    );

    # previous iterate
    X = [deepcopy(x0) for i = 1:N]
    U = [mass * gravity for k = 1:N-1]

    # new iterate
    Xn = deepcopy(X)
    Un = deepcopy(U)

    P = [zeros(nx,nx) for i = 1:N]   # cost to go quadratic term
    p = [zeros(nx) for i = 1:N]      # cost to go linear term
    d = [zeros(nu) for i = 1:N-1]    # feedforward control
    K = [zeros(nu,nx) for i = 1:N-1] # feedback gain
    Xhist = solve!(params,X,U,P,p,K,d,Xn,Un;atol=1e-1,max_iters = 15,verbose = true,ρ = 1e0, ϕ = 10.0);
    function mat_from_vec(X::Vector{Vector{Float64}})::Matrix
        # convert a vector of vectors to a matrix 
        Xm = hcat(X...)
        return Xm 
    end
    Xsim_m = mat_from_vec(Xn)
    Usim_m = mat_from_vec(Un)
    using Plots
    display(plot(t_vec,Xsim_m',label = ["x₁" "x₂" "x₃" "ẋ₁" "ẋ₂" "ẋ₃"],
                title = "State History",
                xlabel = "time (s)", ylabel = "x"))
    display(plot(t_vec[1:end-1],Usim_m',label = ["u₁" "u₂" "u₃"],
                title = "Input History",
                xlabel = "time (s)", ylabel = "u"))
end

main()