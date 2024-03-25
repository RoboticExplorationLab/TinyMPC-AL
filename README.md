# TinyMPC-AL
Under Development

This is a previous attempt of TinyMPC that uses augmented Lagrange method and solve TVLQR online (applicable to nonlinear dynamics). No tricks are played here so it is slow.

New developments are at [https://tinympc.org/](https://tinympc.org/)

## Descriptions

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

## How to compile and run (tested on Linux)

1. Clone this repo

```bash
git clone https://github.com/RoboticExplorationLab/TinyMPC-AL.git
```

2. Build the source code

```bash
cd TinyMPC-AL
cmake -S. -Bbuild
cd build
make
```

3. Run the bicycle example

```bash
./examples/bicycle_example
```

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
- Experiments on double integrator, planar quadrotor and bicycle model.
