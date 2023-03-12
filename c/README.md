# README

- AL-TVLQR is ready to use. It is able to handle input/state box constraints and
goal constraint within stabilization or tracking problems for LTI and LTV
systems. Check `test/al_lqr_test` for all tests, experiments and examples. Check
`examples` for MPC experiments.
- AL-iLQR is under development. It is able to handle input/state box constraints
and goal constraint within stabilization or tracking problems for nonlinear
systems.

## How to compile and run

1. Inside `TinyMPC/c` directory, use `cmake -S. -Bbuild`. This creates a new
`build` directory and configuration.  

2. To build the entire project, use `cmake --build build`. Alternatively, you
can build a particular target by `cmake --build build -t target_name`. Let's
try `cmake --build build -t al_lqr_lti_test`.  

3. To run `al_lqr_lti_test` executable, use `./build/test/al_lqr_test/al_lqr_lti_test` while
still inside `TinyMPC/c` directory.  

4. Iterate via 2 to develop your programs.  

## Notes

- Use `slap_MatMulAdd(C, A, B, 1, 0)` instead of `slap_MatMulAB(C, A, B)`
because the later one ignores all metadata.  
- Can use `slap_MatrixAddition(C, C, A, alp)` to bias `C = C + alp*A`.  
- Should use zero-initialization of array in global scope.  
- Should pass by reference instead of return type  
- Linear term q, qf, r come from reference trajectories, ie. q = -Q*xref

## Done

- Augmented Lagrange LQR for LTV systems.
- Test all units and integration.
- Experiment with double integrator, planar quadrotor and bicycle model.
