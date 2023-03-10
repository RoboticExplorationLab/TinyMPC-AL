#include "constrained_lqr.h"

void tiny_AddStageCost(double* cost, const tiny_ProblemData prob,
                       const Matrix x, const Matrix u, const int k) {
  double dx_data[prob.nstates];
  Matrix dx = slap_MatrixFromArray(prob.nstates, 1, dx_data);
  slap_MatrixAddition(dx, x, prob.X_ref[k], -1);
  *cost += 0.5 * slap_QuadraticForm(dx, prob.Q, dx);
  Matrix du = slap_MatrixFromArray(prob.ninputs, 1, dx_data);
  slap_MatrixAddition(du, u, prob.U_ref[k], -1);
  *cost += 0.5 * slap_QuadraticForm(du, prob.R, du);
}

void tiny_AddTerminalCost(double* cost, const tiny_ProblemData prob,
                          const Matrix x) {
  double dx_data[prob.nstates];
  Matrix dx = slap_MatrixFromArray(prob.nstates, 1, dx_data);
  slap_MatrixAddition(dx, x, prob.X_ref[prob.nhorizon - 1], -1);
  *cost += 0.5 * slap_QuadraticForm(dx, prob.Qf, dx);
}

void tiny_ExpandStageCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                          Matrix* hes_el_uu, Matrix* grad_el_u,
                          const tiny_ProblemData prob, const int k) {
  slap_MatrixCopy(*hes_el_xx, prob.Q);
  slap_MatMulAdd(*grad_el_x, prob.Q, prob.X_ref[k], -1, 0);
  slap_MatrixCopy(*hes_el_uu, prob.R);
  slap_MatMulAdd(*grad_el_u, prob.R, prob.U_ref[k], -1, 0);
}

void tiny_ExpandTerminalCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                             const tiny_ProblemData prob) {
  slap_MatrixCopy(*hes_el_xx, prob.Qf);
  slap_MatMulAdd(*grad_el_x, prob.Qf, prob.X_ref[prob.nhorizon-1], -1, 0);
}

// Riccati recursion for LTI without constraints
enum slap_ErrorCode tiny_BackwardPassLti(tiny_ProblemData* prob,
                                         const tiny_Solver solver,
                                         const tiny_LinearDiscreteModel model,
                                         const Matrix* X, const Matrix* U,
                                         Matrix Q_temp) {
  int N = prob->nhorizon;
  int n = prob->nstates;
  int m = prob->ninputs;
  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob);

  Matrix Qxx = slap_CreateSubMatrix(Q_temp, 0, 0, n, n);
  Matrix Qxu = slap_CreateSubMatrix(Q_temp, 0, n, n, m);
  Matrix Qux = slap_CreateSubMatrix(Q_temp, n, 0, m, n);
  Matrix Quu = slap_CreateSubMatrix(Q_temp, n, n, m, m);
  Matrix Qx = slap_CreateSubMatrix(Q_temp, 0, n + m, n, 1);
  Matrix Qu = slap_CreateSubMatrix(Q_temp, n, n + m, m, 1);

  // Hijack the first part of P[N] to use for Cholesky decomposition
  // NOTE: Assumes m <= n
  Matrix Quu_temp = slap_Reshape(prob->P[N - 1], m, m);
  for (int k = N - 2; k >= 0; --k) {
    // Stage cost expansion
    tiny_ExpandStageCost(&Qxx, &Qx, &Quu, &Qu, *prob, k);
    // State Gradient: Qx = q + A'(P*f + p)
    slap_MatrixCopy(prob->p[k], prob->p[k + 1]);
    // slap_MatMulAdd(prob->p[k], prob->P[k + 1], model.f, 1,
                  //  1);  // p[k] = P[k+1]*f + p[k+1]
    slap_MatMulAdd(Qx, slap_Transpose(model.A), prob->p[k], 1, 1);

    // Control Gradient: Qu = r + B'(P*f + p)
    slap_MatMulAdd(Qu, slap_Transpose(model.B), prob->p[k], 1, 1);

    // State Hessian: Qxx = Q + A'P*A
    slap_MatMulAdd(prob->P[k], prob->P[k + 1], model.A, 1,
                   0);  // P[k] = P[k+1]*A
    slap_MatMulAdd(Qxx, slap_Transpose(model.A), prob->P[k], 1,
                   1);  // Qxx = Q + A'P*A

    // Control Hessian Quu = R + B'P*B
    slap_MatMulAdd(Qxu, prob->P[k + 1], model.B, 1, 0);       // Qxu = P * B
    slap_MatMulAdd(Quu, slap_Transpose(model.B), Qxu, 1, 1);  // Quu = R + B'P*B
    slap_MatMulAdd(Quu, slap_Transpose(model.B), model.B, solver.regu, 1);
    // Hessian Cross-Term
    slap_MatMulAdd(Qux, slap_Transpose(model.B), prob->P[k], 1,
                   0);  // Qux = B'P*A

    // Calculate Gains
    slap_MatrixCopy(Quu_temp, Quu);
    slap_Cholesky(Quu_temp);
    slap_MatrixCopy(prob->K[k], Qux);
    slap_MatrixCopy(prob->d[k], Qu);
    slap_CholeskySolve(Quu_temp, prob->d[k]);  // d = Quu\Qu
    slap_CholeskySolve(Quu_temp, prob->K[k]);  // K = Quu\Qux

    // Cost-to-Go Hessian: P = Qxx + K'Quu*K - K'Qux - Qux'K
    slap_MatrixCopy(prob->P[k], Qxx);                            // P = Qxx
    slap_MatMulAdd(Qxu, slap_Transpose(prob->K[k]), Quu, 1, 0);  // Qxu = K'Quu
    slap_MatMulAdd(prob->P[k], Qxu, prob->K[k], 1, 1);           // P += K'Quu*K
    slap_MatMulAdd(prob->P[k], slap_Transpose(prob->K[k]), Qux, -2,
                   1);  // P -= K'Qux
    // slap_MatMulAdd(prob->P[k], slap_Transpose(Qux), prob->K[k], -1,
    //                1);  // P -= Qux'K

    // Cost-to-Go Gradient: p = Qx + K'Quu*d - K'Qu - Qux'd
    slap_MatrixCopy(prob->p[k], Qx);                    // p = Qx
    slap_MatMulAdd(prob->p[k], Qxu, prob->d[k], 1, 1);  // p += K'Quu*d
    slap_MatMulAdd(prob->p[k], slap_Transpose(prob->K[k]), Qu, -1,
                   1);  // p -= K'Qu
    slap_MatMulAdd(prob->p[k], slap_Transpose(Qux), prob->d[k], -1,
                   1);  // p -= Qux'd
  }
  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob);  // Replace P[N] since we used it for
  // Quu_temp (need improving later)
  return SLAP_NO_ERROR;
}

// // Riccati recursion for LTI with constraints
// enum slap_ErrorCode tiny_ConstrainedBackwardPassLti(
//     tiny_ProblemData* prob, const tiny_Solver solver, 
//     const tiny_LinearDiscreteModel model, const Matrix* X, const Matrix* U, 
//     Matrix* Q_temp, Matrix* ineq_temp) {
//   // Copy terminal cost-to-go
//   int N = prob->nhorizon;
//   tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob,
//                           X[N - 1]);
//   int n = prob->nstates;
//   int m = prob->ninputs;

//   Matrix Qxx = slap_CreateSubMatrix(*Q_temp, 0, 0, n, n);
//   Matrix Qxu = slap_CreateSubMatrix(*Q_temp, 0, n, n, m);
//   Matrix Qux = slap_CreateSubMatrix(*Q_temp, n, 0, m, n);
//   Matrix Quu = slap_CreateSubMatrix(*Q_temp, n, n, m, m);
//   Matrix Qx = slap_CreateSubMatrix(*Q_temp, 0, n + m, n, 1);
//   Matrix Qu = slap_CreateSubMatrix(*Q_temp, n, n + m, m, 1);

//   // Hijack the first part of P[N] to use for Cholesky decomposition
//   // NOTE: Assumes m <= n
//   Matrix Quu_temp = slap_Reshape(prob->P[N - 1], m, m);

//   Matrix ineq_input =
//       slap_CreateSubMatrix(*ineq_temp, 0, 0, prob->ncstr_inputs, 1);
//   Matrix ineq_input2 =
//       slap_CreateSubMatrix(*ineq_temp, 0, 1, prob->ncstr_inputs, 1);
//   Matrix mask_input = slap_CreateSubMatrix(*ineq_temp, 0, 2, prob->ncstr_inputs,
//                                            prob->ncstr_inputs);
//   Matrix ineq_input_jac = slap_CreateSubMatrix(
//       *ineq_temp, 0, prob->ncstr_inputs + 2, prob->ncstr_inputs, prob->ninputs);
//   Matrix ineq_input_jac2 = slap_CreateSubMatrix(
//       *ineq_temp, 0, prob->ncstr_inputs + 2 + prob->ninputs, prob->ncstr_inputs,
//       prob->ninputs);

//   Matrix ineq_state =
//       slap_CreateSubMatrix(*ineq_temp, 0, 0, prob->ncstr_states, 1);
//   Matrix ineq_state2 =
//       slap_CreateSubMatrix(*ineq_temp, 0, 1, prob->ncstr_states, 1);
//   Matrix mask_state = slap_CreateSubMatrix(*ineq_temp, 0, 2, prob->ncstr_states,
//                                            prob->ncstr_states);
//   Matrix ineq_state_jac = slap_CreateSubMatrix(
//       *ineq_temp, 0, prob->ncstr_states + 2, prob->ncstr_states, prob->nstates);
//   Matrix ineq_state_jac2 = slap_CreateSubMatrix(
//       *ineq_temp, 0, prob->ncstr_states + 2 + prob->nstates, prob->ncstr_states,
//       prob->nstates);

//   //========= Goal constraints ==========
//   slap_MatrixAddition(Qx, X[N - 1], prob->X_ref[N - 1], -1);  // hxv
//   slap_SetIdentity(Qxx, 1);                                   // grad_hx
//   slap_MatrixAddition(prob->p[N - 2], prob->goal_dual, Qx,
//                       solver.penalty);  // p[N]  += ∇hx'*(λ + ρ*hxv)
//   slap_MatMulAdd(prob->p[N - 1], slap_Transpose(Qxx), prob->p[N - 2], 1,
//                  1);  // !!constant jacobian
//   slap_MatMulAdd(prob->P[N - 1], slap_Transpose(Qxx), Qxx, solver.penalty,
//                  1);  // P[N]  += ρ*∇hx'∇hx

//   //========= State constraints at end ==========
//   tiny_IneqStates(&ineq_state, *prob, X[N - 1]);  // ineq_state size = 2*NINPUTS
//   tiny_ActiveIneqMask(&mask_state, prob->state_duals[N - 1], ineq_state);
//   slap_ScaleByConst(mask_state, solver.penalty);  // mask = ρ*mask
//   tiny_IneqStatesJacobian(&ineq_state_jac, *prob);
//   // Qx  += ∇hx'*(μx[k] + ρ*(mask * hxv))
//   slap_MatrixCopy(ineq_state2, prob->state_duals[N - 1]);
//   slap_MatMulAdd(ineq_state2, mask_state, ineq_state, 1, 1);
//   slap_MatMulAdd(prob->p[N - 1], slap_Transpose(ineq_state_jac), ineq_state2, 1,
//                  1);
//   // Qxx += ρ*∇hx'*mask*∇hx
//   slap_MatMulAdd(ineq_state_jac2, mask_state, ineq_state_jac, 1, 0);
//   slap_MatMulAdd(prob->P[N - 1], slap_Transpose(ineq_state_jac),
//                  ineq_state_jac2, 1, 1);

//   for (int k = N - 2; k >= 0; --k) {
//     // Stage cost expansion
//     tiny_ExpandStageCost(&Qxx, &Qx, &Quu, &Qu, *prob, X[k], U[k], k);
//     // State Gradient: Qx = q + A'(P*f + p)
//     slap_MatrixCopy(prob->p[k], prob->p[k + 1]);
//     // slap_MatMulAdd(prob->p[k], prob->P[k + 1], model.f, 1,
//     //                1);  // p[k] = P[k+1]*f + p[k+1]
//     slap_MatMulAdd(Qx, slap_Transpose(model.A), prob->p[k], 1, 1);

//     // Control Gradient: Qu = r + B'(P*f ineq_input_jac+ p)
//     slap_MatMulAdd(Qu, slap_Transpose(model.B), prob->p[k], 1, 1);

//     // State Hessian: Qxx = Q + A'P*A
//     slap_MatMulAdd(prob->P[k], prob->P[k + 1], model.A, 1,
//                    0);  // P[k] = P[k+1]*A
//     slap_MatMulAdd(Qxx, slap_Transpose(model.A), prob->P[k], 1,
//                    1);  // Qxx = Q + A'P*A

//     // Control Hessian Quu = R + B'P*B
//     slap_MatMulAdd(Qxu, prob->P[k + 1], model.B, 1, 0);       // Qxu = P * B
//     slap_MatMulAdd(Quu, slap_Transpose(model.B), Qxu, 1, 1);  // Quu = R + B'P*B
//     slap_MatMulAdd(Quu, slap_Transpose(model.B), model.B, solver.regu, 1);
//     // Hessian Cross-Term
//     slap_MatMulAdd(Qux, slap_Transpose(model.B), prob->P[k], 1,
//                    0);  // Qux = B'P*A

//     //========= Control constraints ==========
//     tiny_IneqInputs(&ineq_input, *prob, U[k]);  // ineq_input size = 2*NINPUTS
//     tiny_ActiveIneqMask(&mask_input, prob->input_duals[k], ineq_input);
//     slap_ScaleByConst(mask_input, solver.penalty);  // mask = ρ*mask
//     tiny_IneqInputsJacobian(&ineq_input_jac, *prob);
//     // Qu  += ∇hu'*(μ[k] + (mask * huv))
//     slap_MatrixCopy(ineq_input2, prob->input_duals[k]);
//     slap_MatMulAdd(ineq_input2, mask_input, ineq_input, 1, 1);
//     slap_MatMulAdd(Qu, slap_Transpose(ineq_input_jac), ineq_input2, 1, 1);
//     // Quu += ∇hu'*mask*∇hu
//     slap_MatMulAdd(ineq_input_jac2, mask_input, ineq_input_jac, 1, 0);
//     slap_MatMulAdd(Quu, slap_Transpose(ineq_input_jac), ineq_input_jac2, 1, 1);

//     //========= State constraints ==========
//     tiny_IneqStates(&ineq_state, *prob, X[k]);  // ineq_state size = 2*NINPUTS
//     tiny_ActiveIneqMask(&mask_state, prob->state_duals[k], ineq_state);
//     slap_ScaleByConst(mask_state, solver.penalty);  // mask = ρ*mask
//     tiny_IneqStatesJacobian(&ineq_state_jac, *prob);
//     // Qx  += ∇hx'*(μx[k] + ρ*(mask * hxv))
//     slap_MatrixCopy(ineq_state2, prob->state_duals[k]);
//     slap_MatMulAdd(ineq_state2, mask_state, ineq_state, 1, 1);
//     slap_MatMulAdd(Qx, slap_Transpose(ineq_state_jac), ineq_state2, 1, 1);
//     // Qxx += ρ*∇hx'*mask*∇hx
//     slap_MatMulAdd(ineq_state_jac2, mask_state, ineq_state_jac, 1, 0);
//     slap_MatMulAdd(Qxx, slap_Transpose(ineq_state_jac), ineq_state_jac2, 1, 1);

//     // Calculate Gains
//     slap_MatrixCopy(Quu_temp, Quu);
//     slap_Cholesky(Quu_temp);
//     slap_MatrixCopy(prob->K[k], Qux);
//     slap_MatrixCopy(prob->d[k], Qu);
//     slap_CholeskySolve(Quu_temp, prob->d[k]);  // d = Quu\Qu
//     slap_CholeskySolve(Quu_temp, prob->K[k]);  // K = Quu\Qux

//     // Cost-to-Go Hessian: P = Qxx + K'Quu*K - K'Qux - Qux'K
//     slap_MatrixCopy(prob->P[k], Qxx);                            // P = Qxx
//     slap_MatMulAdd(Qxu, slap_Transpose(prob->K[k]), Quu, 1, 0);  // Qxu = K'Quu
//     slap_MatMulAdd(prob->P[k], Qxu, prob->K[k], 1, 1);           // P += K'Quu*K
//     slap_MatMulAdd(prob->P[k], slap_Transpose(prob->K[k]), Qux, -2,
//                    1);  // P -= 2*K'Qux
//     // slap_MatMulAdd(prob->P[k], slap_Transpose(Qux), prob->K[k], -1,
//     //                1);  // P -= Qux'K

//     // Cost-to-Go Gradient: p = Qx + K'Quu*d - K'Qu - Qux'd
//     slap_MatrixCopy(prob->p[k], Qx);                    // p = Qx
//     slap_MatMulAdd(prob->p[k], Qxu, prob->d[k], 1, 1);  // p += K'Quu*d
//     slap_MatMulAdd(prob->p[k], slap_Transpose(prob->K[k]), Qu, -1,
//                    1);  // p -= K'Qu
//     slap_MatMulAdd(prob->p[k], slap_Transpose(Qux), prob->d[k], -1,
//                    1);  // p -= Qux'd
//   }

//   tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob,
//                           X[N - 1]);  // Replace P[N] since we used it for
//                                       // Quu_temp (need improving later)
//   return SLAP_NO_ERROR;
// }

// enum slap_ErrorCode tiny_AugmentedLagrangianLqr(
//     Matrix* X, Matrix* U, tiny_ProblemData* prob, tiny_Solver* solver,
//     const tiny_LinearDiscreteModel model, const int verbose) {
//   int N = prob->nhorizon;
//   int n = prob->nstates;
//   int m = prob->ninputs;
//   for (int k = 0; k < N - 1; ++k) {
//     tiny_DiscreteDynamics(&(X[k + 1]), X[k], U[k], model);
//   }
//   double G_temp_data[(n + m) * (n + m + 1)];
//   Matrix Q_temp = slap_MatrixFromArray(n + m, n + m + 1, G_temp_data);

//   double ineq_temp_data[prob->ncstr_states *
//                         (prob->ncstr_states + prob->ncstr_states + 2)];
//   Matrix ineq_temp = slap_MatrixFromArray(
//       prob->ncstr_states, prob->ncstr_states + 2 * prob->nstates + 2,
//       ineq_temp_data);

//   Matrix ineq_input =
//       slap_CreateSubMatrix(ineq_temp, 0, 0, prob->ncstr_inputs, 1);
//   Matrix new_input_duals =
//       slap_CreateSubMatrix(ineq_temp, 0, 1, prob->ncstr_inputs, 1);
//   Matrix mask_input = slap_CreateSubMatrix(ineq_temp, 0, 2, prob->ncstr_inputs,
//                                            prob->ncstr_inputs);

//   Matrix ineq_state =
//       slap_CreateSubMatrix(ineq_temp, 0, 0, prob->ncstr_states, 1);
//   Matrix new_state_duals =
//       slap_CreateSubMatrix(ineq_temp, 0, 1, prob->ncstr_states, 1);
//   Matrix mask_state = slap_CreateSubMatrix(ineq_temp, 0, 2, prob->ncstr_states,
//                                            prob->ncstr_states);

//   Matrix eq_goal = slap_MatrixFromArray(prob->ncstr_goal, 1, ineq_temp_data);

//   double norm_d_max = 0.0;
//   double cstr_violation = 0.0;
//   for (int iter = 0; iter < solver->max_primal_iters; ++iter) {
//     tiny_ConstrainedBackwardPassLti(prob, model, *solver, X, U, &Q_temp,
//                                     &ineq_temp);
//     tiny_ForwardPassLti(X, U, *prob, model);

//     norm_d_max = tiny_RiccatiConvergence(*prob);

//     if (verbose == 1) {
//       printf("outer loop\n");
//       printf(
//           "iter     J           ΔJ        |d|         α        reg         "
//           "ρ\n");
//       printf(
//           "--------------------------------------------------------------------"
//           "-\n");
//       printf("%3d   %10.3e  %9.2e  %9.2e  %6.4f   %9.2e   %9.2e\n", iter, 0.0,
//              0.0, norm_d_max, 1.0, solver->regu, solver->penalty);
//     }

//     printf("update duals and penalty\n");

//     // For linear systems, only 1 iteration, shouldn't need condition here
//     if (0*norm_d_max < solver->riccati_tol) {
//       cstr_violation = 0.0;
//       double norm_inf = 0.0;

//       for (int k = 0; k < N - 1; ++k) {
//         //========= Control constraints ==========
//         tiny_IneqInputs(&ineq_input, *prob,
//                         U[k]);  // ineq_input size = 2*NINPUTS
//         tiny_ActiveIneqMask(&mask_input, prob->input_duals[k], ineq_input);
//         slap_ScaleByConst(mask_input, solver->penalty);  // mask = ρ*mask
//         // Update duals
//         slap_MatrixCopy(new_input_duals, prob->input_duals[k]);
//         slap_MatMulAdd(new_input_duals, mask_input, ineq_input, 1,
//                        1);  //μ[k] + ρ*mask*huv
//         slap_ArgMax(ineq_input, &norm_inf);
//         norm_inf = norm_inf > 0.0 ? norm_inf : 0.0;
//         norm_inf = norm_inf * 2;
//         tiny_ClampIneqDuals(&(prob->input_duals)[k], new_input_duals);
//         // convio = max(convio,norm(hxv + abs.(hxv),Inf))
//         cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
//       }

//       for (int k = 0; k < N; ++k) {
//         //========= State constraints ==========
//         tiny_IneqStates(&ineq_state, *prob,
//                         X[k]);  // ineq_input size = 2*NINPUTS
//         tiny_ActiveIneqMask(&mask_state, prob->state_duals[k], ineq_state);
//         slap_ScaleByConst(mask_state, solver->penalty);  // mask = ρ*mask
//         // Update duals
//         slap_MatrixCopy(new_state_duals, prob->state_duals[k]);
//         slap_MatMulAdd(new_state_duals, mask_state, ineq_state, 1,
//                        1);  //μ[k] + ρ*mask*huv
//         slap_ArgMax(ineq_state, &norm_inf);
//         norm_inf = norm_inf > 0.0 ? norm_inf : 0.0;
//         norm_inf = norm_inf * 2;
//         tiny_ClampIneqDuals(&(prob->state_duals)[k], new_state_duals);
//         // convio = max(convio,norm(hxv + abs.(hxv),Inf))
//         cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
//       }

//       //========= Goal constraints ==========
//       slap_MatrixAddition(eq_goal, X[N - 1], prob->X_ref[N - 1], -1);
//       norm_inf = slap_NormInf(eq_goal);
//       cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
//       // λ .+= ρ*hxv
//       slap_MatrixAddition(prob->goal_dual, prob->goal_dual, eq_goal,
//                           solver->penalty);

//       printf("convio: %.6f \n", cstr_violation);
//       if (cstr_violation < solver->cstr_tol) {
//         printf("\nSUCCESS!\n");
//         return SLAP_NO_ERROR;
//       }
//       solver->penalty = solver->penalty * solver->penalty_mul;
//     }
//   }
//   return SLAP_NO_ERROR;
// }

// Roll out the closed-loop dynamics with K and d from backward pass,
// calculate new X, U in place
enum slap_ErrorCode tiny_ForwardPassLti(Matrix* X, Matrix* U,
                                        const tiny_ProblemData prob,
                                        const tiny_LinearDiscreteModel model) {
  int N = prob.nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    // delta_x and delta_u over previous X, U
    // Control input: u = uf - d - K*(x - xf)
    slap_MatrixCopy(U[k], prob.d[k]);              // u[k] = d[k]
    slap_MatMulAdd(U[k], prob.K[k], X[k], -1, -1);   // u[k] += K[k] * x[k]
    // Next state: x = A*x + B*u + f
    tiny_DiscreteDynamics(&X[k + 1], X[k], U[k], model);
  }
  return SLAP_NO_ERROR;
}

double tiny_RiccatiConvergence(const tiny_ProblemData prob) {
  double norm_d_max = 0.0;
  for (int k = 0; k < prob.nhorizon - 1; ++k) {
    double norm_d = slap_NormTwo(prob.d[k]);
    if (norm_d > norm_d_max) {
      norm_d_max = norm_d;
    }
  }
  return norm_d_max;
}

// [u-p.u_max;-u + p.u_min]
void tiny_IneqInputs(Matrix* ineq_input, const tiny_ProblemData prob,
                     const Matrix u) {
  slap_SetConst(*ineq_input, 0);  // clear before processing
  Matrix upper_half = slap_CreateSubMatrix(*ineq_input, 0, 0, prob.ninputs, 1);
  Matrix lower_half =
      slap_CreateSubMatrix(*ineq_input, prob.ninputs, 0, prob.ninputs, 1);
  slap_MatrixAddition(upper_half, u, prob.u_max, -1);
  slap_MatrixAddition(lower_half, prob.u_min, u, -1);
}

// [u_max, -u_min]
void tiny_IneqInputsOffset(Matrix* ineq_input, const tiny_ProblemData prob) {
  Matrix upper_half = slap_CreateSubMatrix(*ineq_input, 0, 0, prob.ninputs, 1);
  Matrix lower_half =
      slap_CreateSubMatrix(*ineq_input, prob.ninputs, 0, prob.ninputs, 1);
  slap_MatrixCopy(upper_half, prob.u_max);
  slap_MatrixCopy(lower_half, prob.u_min);
  slap_ScaleByConst(lower_half, -1);
}

void tiny_IneqInputsJacobian(Matrix* ineq_jac, const tiny_ProblemData prob) {
  slap_SetConst(*ineq_jac, 0);  // clear before processing
  Matrix upper_half =
      slap_CreateSubMatrix(*ineq_jac, 0, 0, prob.ninputs, prob.ninputs);
  Matrix lower_half = slap_CreateSubMatrix(*ineq_jac, prob.ninputs, 0,
                                           prob.ninputs, prob.ninputs);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
}

// [x-p.x_max;-x + p.x_min]
void tiny_IneqStates(Matrix* ineq_state, const tiny_ProblemData prob,
                     const Matrix x) {
  slap_SetConst(*ineq_state, 0);  // clear before processing
  Matrix upper_half = slap_CreateSubMatrix(*ineq_state, 0, 0, prob.nstates, 1);
  Matrix lower_half =
      slap_CreateSubMatrix(*ineq_state, prob.nstates, 0, prob.nstates, 1);
  slap_MatrixAddition(upper_half, x, prob.x_max, -1);
  slap_MatrixAddition(lower_half, prob.x_min, x, -1);
}

// [x_max, -x_min]
void tiny_IneqStatesOffset(Matrix* ineq_state, const tiny_ProblemData prob) {
  Matrix upper_half = slap_CreateSubMatrix(*ineq_state, 0, 0, prob.nstates, 1);
  Matrix lower_half =
      slap_CreateSubMatrix(*ineq_state, prob.nstates, 0, prob.nstates, 1);
  slap_MatrixCopy(upper_half, prob.x_max);
  slap_MatrixCopy(lower_half, prob.x_min);
  slap_ScaleByConst(lower_half, -1);
} 

void tiny_IneqStatesJacobian(Matrix* ineq_jac, const tiny_ProblemData prob) {
  slap_SetConst(*ineq_jac, 0);  // clear before processing
  Matrix upper_half =
      slap_CreateSubMatrix(*ineq_jac, 0, 0, prob.nstates, prob.nstates);
  Matrix lower_half = slap_CreateSubMatrix(*ineq_jac, prob.nstates, 0,
                                           prob.nstates, prob.nstates);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
}


void tiny_ActiveIneqMask(Matrix* mask, const Matrix dual, const Matrix ineq) {
  slap_SetConst(*mask, 0);  // clear before processing
  for (int i = 0; i < ineq.rows; ++i) {
    // When variables are on the boundary or violating constraints
    bool active = dual.data[i] > 0 || ineq.data[i] > 0;
    slap_SetElement(*mask, i, i, active);
  }
}

void tiny_ClampIneqDuals(Matrix* dual, const Matrix new_dual) {
  for (int i = 0; i < dual->rows; ++i) {
    if (new_dual.data[i] > 0) {
      dual->data[i] = new_dual.data[i];
    } else
      dual->data[i] = 0.0;
  }
}

void tiny_DiscreteDynamics(Matrix* xn, const Matrix x, const Matrix u,
                           const tiny_LinearDiscreteModel model) {
  slap_MatrixCopy(*xn, model.f);                          
  slap_MatMulAdd(*xn, model.A, x, 1, 1);      // x[k+1] += A * x[k]
  slap_MatMulAdd(*xn, model.B, u, 1, 1);      // x[k+1] += B * u[k]
}