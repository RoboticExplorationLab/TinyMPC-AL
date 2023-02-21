using TrajectoryOptimization
import RobotZoo.Rover
using StaticArrays, LinearAlgebra
using Plots
using RobotDynamics
using Altro
const RD = RobotDynamics

print("LQR SIMULATION START!!!")
# QP, Optimal control note
# Equivalent to LQR

# SET UP VISUALIZATION, SHOULD NOT RERUN
# Create Visualizer
# vis = rover_vis_init()

# SOLVER
dt = 1/20   # 10 - 20Hz
N = 151
tf = (N-1)*dt

steer_rate_lim=0.3
accel_lim=3

model = Rover(ref=:rear)
n,m = size(model);

x0 = @SVector zeros(n)
xf = @SVector [15, 15, π/4, 0, 0];  

# Objective function
Q = 1e-1*Diagonal(@SVector ones(n)) 
Qf = 1e3*Diagonal(@SVector ones(n))
R = 1e1*Diagonal(@SVector ones(m)) 
obj = LQRObjective(Q,R,Qf,xf,N);

# Create Empty ConstraintList
conSet = ConstraintList(n,m,N)

# Circle constraints
xc = SA[6, 10, 12]
yc = SA[5, 12, 5]
rc = SA[2, 2, 2]    #try 1
cir = CircleConstraint(n, xc, yc, rc)
add_constraint!(conSet, cir, 2:N-1)

# Control-state Bounds
u_bnd = [accel_lim; steer_rate_lim]
x_bnd = [Inf;Inf;Inf;40;π/3]
bnd = BoundConstraint(n,m, x_min=-x_bnd, x_max=x_bnd, u_min=-u_bnd, u_max=u_bnd)
add_constraint!(conSet, bnd, 1:N-1)
# Goal Constraint
goal = GoalConstraint(xf,1:2)
add_constraint!(conSet, goal, N)

prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)

u0 = @SVector fill(0.0,m)
U0 = [u0 for k = 1:N-1]
initial_controls!(prob, U0)
rollout!(prob);

solver = ALTROSolver(prob)
set_options!(solver, show_summary=true)
initial_controls!(solver, U0)
solve!(solver);

Usln = hcat(Vector.(controls(solver))...)
Xsln = hcat(Vector.(states(solver))...)

# Rover simulator
xs = zeros(5,N)
x0_off = [0;2;0;0;0]
xs[:,1] .= x0 + x0_off
dd = zeros(5,N)

# dmodel2 = RD.DiscretizedDynamics{RD.RK4}(model2)
for k in 1:N-1
    dd[:,k] = 0*randn(5)
    # xs[:,k+1] .= simulator(model, xs[:,k]+dd[:,k], Usln[:,k], u_bnd, dt)
    xs[:,k+1] .= DT_non_model(model, xs[:,k]+dd[:,k], Usln[:,k], dt)
end

# Got the reference trajectory, solve for controls with LQR

thist = Array(range(0,dt*(N-1), step=dt))
Un = copy(Usln)  # ref control
Xn = copy(Xsln)  # ref state
Q = Array(1*I(5))
R = Array(1*I(2))
Qf = Array(1.0*I(5))

# Use discrete TV-LQR with Riccati backward recursion

P = zeros(n,n,N)
K = zeros(m,n,N-1)
P[:,:,N] .= Qf
for k = (N-1):-1:1
    AA,BB = DT_lin_model(model, Xn[:,k], Un[:,k], dt)
    K[:,:,k] .= (R + BB'*P[:,:,k+1]*BB)\(BB'*P[:,:,k+1]*AA)
    P[:,:,k] .= Q + AA'*P[:,:,k+1]*(AA - BB*K[:,:,k])
end

#Forward rollout starting at x0
xhist = zeros(n,N)
xhist[:,1] = xs[:,1]
uhist = zeros(m,N-1)

for k = 1:(N-1)
    uhist[:,k] .= Un[:,k] - K[:,:,k,1]*(xhist[:,k] - Xn[:,k])
    xhist[:,k+1] .= DT_non_model(model, xhist[:,k]+dd[:,k], uhist[:,k], dt)
    # Ak, Bk = DT_lin_model(model,Xn[:,k],Un[:,k],dt)
    # xhist[:,k+1] .= Xn[:,k+1] + Ak*(xhist[:,k] - Xn[:,k]) + Bk*(uhist[:,k] - Un[:,k])
end


# Plot solution
layout = @layout [a; b]

p1 = plot(xs[1,1:end], xs[2,1:end], xlabel="x position", ylabel="y position", color="gray",linewidth=3,widen=1.1,label="ALTRO")
plot!(Xn[1,1:end], Xn[2,1:end], color="green",linewidth = 5,linestyle=:dot,label="nominal")
scatter!([x0[1]],[x0[2]],color = "yellow", label = "start", markershape=:rect, markersize = 10,alpha = 0.8)
scatter!([xf[1]],[xf[2]],color = "green", label = "goal", markershape=:star, markersize = 12,alpha = 0.8)
scatter!(xc, yc, markershape=:circle, label="obs", color ="red",markersize=rc*15,legend=:bottomright)
plot!(xhist[1,1:end], xhist[2,1:end], color="blue",linewidth=3,label="QP")
# Plot control inputs
time = gettimes(solver)[1:end-1]
p2 = plot(time, Usln[1,1:end], label="accel",lw=4)
plot!(time, uhist[1,1:end], label="accel+",lw=2)
plot!(time, 10*Usln[2,1:end], label="10*steer", xlabel="time (s)", lw=4)
plot!(time, 10*uhist[2,1:end], label="10*steer+", xlabel="time (s)", ylabel="control input", lw=2, reuse=false)

display(plot(p1, p2, layout=layout, size=(500,600)))

# plot(thist, Xsln[2,:]-xhist[2,:])
#
# Animation
# rover_vis_run(vis, xhist, xc, yc, rc, xf)

print("LQR SIMULATION FINISHED!!!")
