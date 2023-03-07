## How to compile and run
1. `cd build`: create and enter a `build` directory  

2. `cmake ..`: config compiler out-source  

3. `make`: build all targets, or  

`make main_tvlqr`: build a particular one (check `CMakeLists.txt`)

4. `./examples/riccati/main_tvlqr`: still inside `build`, run the executable

5. Iterate via 3 to develop your programs.

## Notes:

- Use `slap_MatMulAdd(C, A, B, 1, 0)` instead of `slap_MatMulAB(C, A, B)` because
the later one ignores all metadata.  
- Can use `slap_MatrixAddition(C, C, A, alp)` to bias `C = C + alp*A`.  
- Should use zero-initialization of array.

## Done:

- LQR for LTI systems with arbitrary start and goal (double integrator).  
- Tracking LQR for LTI systems.  
- LQR for LTV systems (Jacobians fixed) with arbitrary start and goal.  
- LQR for LTV systems (Jacobians compute) with arbitrary start and goal (planar
quadrotor).  
- Complete tracking LQR for LTV systems (bicycle). 