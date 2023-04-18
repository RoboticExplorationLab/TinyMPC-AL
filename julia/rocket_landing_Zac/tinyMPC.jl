# -------------------THIS IS ALL ALTRO-------------------------------
# This enables goal equality constraints, state and input ineq constraints.
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
function conic_cost_expansion(p::NamedTuple, u, λc, ρ, k)
    Uc = cone_u(params, u)
    λhat = 1 * λc # already updated
    # λhat = 1*projection(λc - Uc*cone_scale*ρ)
    ∇c = cone_u_jac(p, u)
    G = -∇c' * ∇projection(λhat) * projection(λhat)
    # H = ∇c'*∇projection(λhat)'*∇projection(λhat)*∇c*ρ
    H = ∇c' * (∇projection(λhat)' * ∇projection(λhat) + ∇²projection(λhat, projection(λhat))) * ∇c * ρ
    return (G - H * params.Uref[k]), H
end
function backward_pass!(params, X, U, P, p, d, K, reg, μ, μx, ρ, λ, λc)
    """LQR backward pass with AL. This !function updates its parameters
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
        if params.ncu_cone > 0
            # AL terms for cone
            Uc = cone_u(params, U[k])
            λhat = 1 * projection(λc[k] - Uc * cone_scale * ρ)
            J += (0.5 / ρ) * (λhat' * λhat - λc[k]' * λc[k])
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
function tiny_solve!(params, X, U, P, p, K, d, Xn, Un; atol=1e-3, max_iters=50, max_inner_iters=10, verbose=true, ρ=1, ρ_max=1e8, ϕ=10)

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
function cone_u(p, u)
    [p.A_cone * u; p.c_cone' * u]
end
function cone_u_jac(p, u)
    J = zeros(p.nu, p.nu)
    J[1:end-1, 1:end] .= p.A_cone
    J[end, 1:end] .= p.c_cone
    return J
end
function ineq_con_u_jac(params, u)
    ForwardDiff.jacobian(_u -> ineq_con_u(params, _u), u)  # lazy way
end
function ineq_con_x_jac(p, x)
    ForwardDiff.jacobian(_x -> ineq_con_x(p, _x), x)  # lazy way
end
function projection(x)
    # assumes x is stacked [v; s] such that ||v||₂ ≤ s
    n = length(x)
    v = view(x, 1:n-1)
    s = x[end]
    a = norm(v)
    if a <= -s          # below the cone 
        return zero(x)
    elseif a <= s       # in the cone
        return x
    elseif a >= abs(s)  # outside the cone 
        return 0.5 * (1 + s / a) * [v; a]
    else
        throw(ErrorException("Invalid second-order cone projection"))  # Nan
    end
end
function ∇projection(x)
    n = length(x)
    J = zeros(eltype(x), n, n)
    s = x[end]
    v = view(x, 1:n-1)
    a = norm(v)
    if a <= -s
        return J  # zeros
    elseif a <= s
        J .= I(n)
        return J  # identity
    elseif a >= abs(s)
        c = 0.5 * (1 + s / a)

        # dvdv ok!
        for i = 1:n-1, j = 1:n-1
            J[i, j] = -0.5 * s / a^3 * v[i] * v[j]
            if i == j
                J[i, j] += c
            end
        end

        # dvds ok!
        for i = 1:n-1
            J[i, n] = 0.5 * v[i] / a
        end

        # dsdv ok!
        for i = 1:n-1
            J[n, i] = ((-0.5 * s / a^2) + c / a) * v[i]
        end
        J[n, n] = 0.5  # ok
        return J
    else
        error("Invalid second-order cone projection.")
    end
    return J
end
function ∇projection(x)
    n = length(x)
    v = view(x, 1:n-1)
    s = x[end]
    a = norm(v)
    J = zeros(n, n)
    if a <= -s                               # below cone
        J .*= 0
    elseif a <= s                            # in cone
        J .*= 0
        for i = 1:n
            J[i, i] = 1.0
        end
    elseif a >= abs(s)                       # outside cone
        # scalar
        c = 0.5 * (1 + s / a)

        # dvdv = dbdv * v' + c * oneunit(SMatrix{n-1,n-1,T})
        for i = 1:n-1, j = 1:n-1
            J[i, j] = -0.5 * s / a^3 * v[i] * v[j]
            if i == j
                J[i, j] += c
            end
        end

        # dvds
        for i = 1:n-1
            J[i, n] = 0.5 * v[i] / a
        end

        # ds
        for i = 1:n-1
            J[n, i] = ((-0.5 * s / a^2) + c / a) * v[i]
        end
        J[n, n] = 0.5
    else
        throw(ErrorException("Invalid second-order cone projection"))
    end
    return J
end
function ∇²projection(x, b)
    # x is lamda_bar, b is projection(lambda_bar)
    n = length(x)
    hess = zeros(eltype(x), n, n)
    v = view(x, 1:n-1)
    bv = view(b, 1:n-1)

    # @assert size(hess) == (n+1,n+1)
    s = x[end]
    bs = b[end]
    a = norm(v)
    vbv = dot(v, bv)

    if a <= -s
        return hess .= 0
    elseif a <= s
        return hess .= 0
    elseif a > abs(s)
        # Original equations from chain rule
        # dvdv = -s/norm(v)^2/norm(v)*(I - (v*v')/(v'v))*bv*v' + 
        #     s/norm(v)*((v*(v'bv))/(v'v)^2 * 2v' - (I*(v'bv) + v*bv')/(v'v)) + 
        #     bs/norm(v)*(I - (v*v')/(v'v))
        # dvds = 1/norm(v)*(I - (v*v')/(v'v))*bv;
        # # display(dvds)
        # # display(dvdv)
        # hess[1:n-1,1:n-1] .= dvdv*0.5
        # hess[1:n-1,n] .= dvds*0.5
        # hess[n:n,1:n-1] .= 0.5*dvds'
        # hess[n,n] = 0
        # return hess

        # The following is just an unrolled version of the above
        n = n - 1
        dvdv = view(hess, 1:n, 1:n)
        dvds = view(hess, 1:n, n + 1)
        dsdv = view(hess, n + 1, 1:n)
        @inbounds for i = 1:n
            hi = 0
            @inbounds for j = 1:n
                Hij = -v[i] * v[j] / a^2
                if i == j
                    Hij += 1
                end
                hi += Hij * bv[j]
            end
            dvds[i] = hi / 2a
            dsdv[i] = dvds[i]
            @inbounds for j = 1:i
                vij = v[i] * v[j]
                H1 = hi * v[j] * (-s / a^3)
                H2 = vij * (2 * vbv) / a^4 - v[i] * bv[j] / a^2
                H3 = -vij / a^2
                if i == j
                    H2 -= vbv / a^2
                    H3 += 1
                end
                H2 *= s / a
                H3 *= bs / a
                dvdv[i, j] = (H1 + H2 + H3) / 2
                dvdv[j, i] = dvdv[i, j]
            end
        end
        hess[end, end] = 0
        return hess
    else
        throw(ErrorException("Invalid second-order cone projection"))
    end
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