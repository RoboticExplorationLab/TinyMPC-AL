# NOTES

## How to compile and run

1. Inside `TinyMPC/c` directory, use `cmake -S. -Bbuild`. This creates a new
`build` directory and configuration.  

2. To build the entire project, use `cmake --build build`. Alternatively, you
can build a particular target by `cmake --build build -t target_name`. Let's
try `cmake --build build -t al_lqr_lti_test`.  

3. To run `al_lqr_lti_test` executable, use `./build/test/al_lqr_lti_test` while
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

- Augmented Lagrange LQR for LTI systems.
- Test all units and integration.
- Experiment with double integrator.
