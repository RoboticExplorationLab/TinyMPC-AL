# README

- This is a full library, embedded (optimized) version is under development.
However, it aims at highly modular integration. You can just use part of the
sources at your need. Currently, assertion is not present. Users should be
responsible for this during development.
- AL-TVLQR is ready to use. It is able to handle input/state box constraints and
goal constraint within stabilization or tracking problems for LTI and LTV
systems. Check `test/al_lqr_test` for all tests, experiments and examples. Check
`examples` for MPC experiments.
- You can set constraints on/off in `CMakeLists.txt`, more are under development.
- AL-iLQR is under development. It is able to handle input/state box constraints
and goal constraint within stabilization or tracking problems for nonlinear
systems.

## How to compile and run

1. Inside `TinyMPC/c` directory, use `cmake -S. -Bbuild`. This creates a new
`build` directory and configuration.  

2. To build the entire project, use `cmake --build build`. Alternatively, you
can build a particular target by `cmake --build build -t target_name`. Let's
try `cmake --build build -t bicycle_example`.  

3. To run `bicycle_example` executable, use `./build/examples/bicycle_example` while
still inside `TinyMPC/c` directory.  

4. Iterate via 2 to develop your programs.  

## Notes

- Use `slap_MatMulAdd(C, A, B, 1, 0)` instead of `slap_MatMulAB(C, A, B)`
because the later one ignores all metadata.  
- Can use `slap_MatrixAddition(C, C, A, alp)` to bias `C = C + alp*A`.  
- Should use zero-initialization of array in global scope.  
- Should pass by reference instead of return type  
- Linear term q, qf, r come from reference trajectories, ie. q = -Q*xref
- MPC for LTI systems can handel all provided types of constraints.
- Tracking MPC for LTV systems may not handle all due to the strictness.

## Done

- Augmented Lagrangian LQR/TVLQR and MPC.
- Successful unit and integration testing.
- Experiments on sfloat integrator, planar quadrotor and bicycle model.

## Optimal Control Problem

- Check [the report](tinyMPC_Report.pdf) for full derivation.

- Check `examples/bicycle_example.c` to see how to use the library (not so good
API yet).

- There are almost no local variables created inside a function (all defined by user before solve).

- Temporary data should not be modified unless you know exactly what you want to do.

- If you have a LTV system, you'd better derive the formulas to compute A(t), B(t), f(t). We use `Symbolics.jl` to derive it analytically and convert to C code.

*Workflow:*

1. Define all necessary array and convert them to `Matrix`.

2. Define necessary struct: `model`, `problem_data`, `solver`.

3. Assign data to struct.

4. Solve the problem.

*Constraints:*

- Inequality constraints will be identically applied to all timestep, in the form: $Ax \leq b$. For example, input bound:
$$A = [I; -I], \quad b = [umax; -umin]$$

- Equality constraint is goal constraint.

- Enable/disable constraints by setting, e.g.,  `prob.ncstr_inputs` to `1` or `0`
