using SparseArrays, OSQP

# Utility function
speye(N) = spdiagm(ones(N))

# Discrete time model of a quadcopter
Ad = [1.0 0.0 0.0 0.05 0.0 0.0; 
    0.0 1.0 0.0 0.0 0.05 0.0; 
    0.0 0.0 1.0 0.0 0.0 0.05; 
    0.0 0.0 0.0 1.0 0.0 0.0; 
    0.0 0.0 0.0 0.0 1.0 0.0; 
    0.0 0.0 0.0 0.0 0.0 1.0] |> sparse
Bd = [0.000125 0.0 0.0; 
0.0 0.000125 0.0;
0.0 0.0 0.000125; 
0.005 0.0 0.0; 
0.0 0.005 0.0; 
0.0 0.0 0.005] |> sparse
(nx, nu) = size(Bd)

# Constraints
u0 = 10.5916
# Sloppy bound to test
umin = -10.0*ones(nu)
umax =  105.0*ones(nu)

# Sloppy bound to test
xmin = [-5,-5,0,-10,-10,-10.0]
xmax = [5,5,20,10,10,10.0]

# Objective function
Q = 10e-1 * speye(nx)
QN = Q
R = 1 * speye(nu)

# Initial and reference states
x0 = [4, 2, 20, -3, 2, -5.0]
xg = [0,0,0,0,0,0.0]
xr = xg
ur = [zeros(nu) for i = 1:N-1] # unused

# Prediction horizon
N = 301

# Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
# - quadratic objective
P = blockdiag(kron(speye(N), Q), QN, kron(speye(N), R))
# - linear objective
q = [repeat(-Q * xr, N); -QN * xr; zeros(N*nu)]
# - linear dynamics
Ax = kron(speye(N + 1), -speye(nx)) + kron(spdiagm(-1 => ones(N)), Ad)
Bu = kron([spzeros(1, N); speye(N)], Bd)
Aeq = [Ax Bu]
leq = [-x0; zeros(N * nx)]
ueq = leq
# - input and state constraints
Aineq = speye((N + 1) * nx + N * nu)
lineq = [repeat(xmin, N + 1); repeat(umin, N)]
uineq = [repeat(xmax, N + 1); repeat(umax, N)]
# - OSQP constraints
A, l, u = [Aeq; Aineq], [leq; lineq], [ueq; uineq]

# Create an OSQP model
m = OSQP.Model()

# Setup workspace
OSQP.setup!(m; P=P, q=q, A=A, l=l, u=u, warm_start=true)

# Simulate in closed loop
nsim = 15;
@time for _ in 1 : nsim
    # Solve
    res = OSQP.solve!(m)

    # Check solver status
    if res.info.status != :Solved
        error("OSQP did not solve the problem!")
    end

    # Apply first control input to the plant
    ctrl = res.x[(N+1)*nx+1:(N+1)*nx+nu]
    global x0 = Ad * x0 + Bd * ctrl

    # Update initial state
    l[1:nx], u[1:nx] = -x0, -x0
    OSQP.update!(m; l=l, u=u)
end