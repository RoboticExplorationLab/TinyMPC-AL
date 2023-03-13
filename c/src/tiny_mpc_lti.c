#include "tiny_mpc_lti.h"

// Riccati recursion for LTI with constraints
enum slap_ErrorCode tiny_ConstrainedBackwardPassLti(
    tiny_ProblemData* prob, const tiny_Solver solver, const tiny_LtiModel model,
    const Matrix* X, const Matrix* U, Matrix* Q_temp, Matrix* ineq_temp) {
  // Copy terminal cost-to-go
  int N = prob->nhorizon;
  int n = prob->nstates;
  int m = prob->ninputs;
  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob);

  Matrix Qxx = slap_CreateSubMatrix(*Q_temp, 0, 0, n, n);
  Matrix Qxu = slap_CreateSubMatrix(*Q_temp, 0, n, n, m);
  Matrix Qux = slap_CreateSubMatrix(*Q_temp, n, 0, m, n);
  Matrix Quu = slap_CreateSubMatrix(*Q_temp, n, n, m, m);
  Matrix Qx = slap_CreateSubMatrix(*Q_temp, 0, n + m, n, 1);
  Matrix Qu = slap_CreateSubMatrix(*Q_temp, n, n + m, m, 1);

  // Hijack the first part of P[N] to use for Cholesky decomposition
  // NOTE: Assumes m <= n
  Matrix Quu_temp = slap_Reshape(prob->P[N - 1], m, m);

  Matrix ineq_input =
      slap_CreateSubMatrix(*ineq_temp, 0, 0, prob->ncstr_inputs, 1);
  Matrix ineq_input2 =
      slap_CreateSubMatrix(*ineq_temp, 0, 1, prob->ncstr_inputs, 1);
  Matrix mask_input = slap_CreateSubMatrix(*ineq_temp, 0, 2, prob->ncstr_inputs,
                                           prob->ncstr_inputs);
  Matrix ineq_input_jac = slap_CreateSubMatrix(
      *ineq_temp, 0, prob->ncstr_inputs + 2, prob->ncstr_inputs, prob->ninputs);
  Matrix ineq_input_jac2 = slap_CreateSubMatrix(
      *ineq_temp, 0, prob->ncstr_inputs + 2 + prob->ninputs, prob->ncstr_inputs,
      prob->ninputs);

  Matrix ineq_state =
      slap_CreateSubMatrix(*ineq_temp, 0, 0, prob->ncstr_states, 1);
  Matrix ineq_state2 =
      slap_CreateSubMatrix(*ineq_temp, 0, 1, prob->ncstr_states, 1);
  Matrix mask_state = slap_CreateSubMatrix(*ineq_temp, 0, 2, prob->ncstr_states,
                                           prob->ncstr_states);
  Matrix ineq_state_jac = slap_CreateSubMatrix(
      *ineq_temp, 0, prob->ncstr_states + 2, prob->ncstr_states, prob->nstates);
  Matrix ineq_state_jac2 = slap_CreateSubMatrix(
      *ineq_temp, 0, prob->ncstr_states + 2 + prob->nstates, prob->ncstr_states,
      prob->nstates);

  //========= Goal constraints ==========
  slap_MatrixCopy(Qx, prob->X_ref[N - 1]);  // h = xg
  slap_SetIdentity(Qxx, 1);                 // H constant
  slap_MatrixAddition(prob->p[N - 2], prob->goal_dual, Qx,
                      -solver.penalty);  // (λ - ρ*h)
  slap_MatMulAdd(prob->p[N - 1], slap_Transpose(Qxx), prob->p[N - 2], 1,
                 1);  // p[N]  += H'*(λ - ρ*h)
  slap_MatMulAdd(prob->P[N - 1], slap_Transpose(Qxx), Qxx, solver.penalty,
                 1);  // P[N]  += ρ*H'H

  //========= State constraints at end ==========
  tiny_IneqStates(&ineq_state, *prob, X[N - 1]);  // ineq_state size = 2*NINPUTS
  tiny_ActiveIneqMask(&mask_state, prob->state_duals[N - 1], ineq_state);
  slap_ScaleByConst(mask_state, solver.penalty);  // mask = ρ*mask
  tiny_IneqStatesJacobian(&ineq_state_jac, *prob);
  // Qx  += G'*(μx[k] - ρ*mask * g)
  tiny_IneqStatesOffset(&ineq_state, *prob);  // g
  slap_MatrixCopy(ineq_state2, prob->state_duals[N - 1]);
  slap_MatMulAdd(ineq_state2, mask_state, ineq_state, -1,
                 1);  //μx[k] - ρ*mask*g
  slap_MatMulAdd(prob->p[N - 1], slap_Transpose(ineq_state_jac), ineq_state2, 1,
                 1);
  // Qxx += G'*ρmask*G
  slap_MatMulAdd(ineq_state_jac2, mask_state, ineq_state_jac, 1, 0);
  slap_MatMulAdd(prob->P[N - 1], slap_Transpose(ineq_state_jac),
                 ineq_state_jac2, 1, 1);

  for (int k = N - 2; k >= 0; --k) {
    // Stage cost expansion
    tiny_ExpandStageCost(&Qxx, &Qx, &Quu, &Qu, *prob, k);
    // State Gradient: Qx = q + A'(P*f + p)
    slap_MatrixCopy(prob->p[k], prob->p[k + 1]);
    slap_MatMulAdd(prob->p[k], prob->P[k + 1], model.f, 1,
                   1);  // p[k] = P[k+1]*f + p[k+1]
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

    //========= Control constraints ==========
    tiny_IneqInputs(&ineq_input, *prob, U[k]);  // ineq_input size = 2*NINPUTS
    tiny_ActiveIneqMask(&mask_input, prob->input_duals[k], ineq_input);
    slap_ScaleByConst(mask_input, solver.penalty);  // mask = ρ*mask
    tiny_IneqInputsJacobian(&ineq_input_jac, *prob);
    // Qu  += G'*(μ[k] + (mask * g)
    tiny_IneqInputsOffset(&ineq_input, *prob);  // g
    slap_MatrixCopy(ineq_input2, prob->input_duals[k]);
    slap_MatMulAdd(ineq_input2, mask_input, ineq_input, -1, 1);
    slap_MatMulAdd(Qu, slap_Transpose(ineq_input_jac), ineq_input2, 1, 1);
    // Quu += ∇hu'*mask*∇hu
    slap_MatMulAdd(ineq_input_jac2, mask_input, ineq_input_jac, 1, 0);
    slap_MatMulAdd(Quu, slap_Transpose(ineq_input_jac), ineq_input_jac2, 1, 1);

    //========= State constraints ==========
    tiny_IneqStates(&ineq_state, *prob, X[k]);  // ineq_state size = 2*NINPUTS
    tiny_ActiveIneqMask(&mask_state, prob->state_duals[k], ineq_state);
    slap_ScaleByConst(mask_state, solver.penalty);  // mask = ρ*mask
    tiny_IneqStatesJacobian(&ineq_state_jac, *prob);
    // Qx  += G'*(μx[k] - ρ*mask * g)
    tiny_IneqStatesOffset(&ineq_state, *prob);  // g
    slap_MatrixCopy(ineq_state2, prob->state_duals[k]);
    slap_MatMulAdd(ineq_state2, mask_state, ineq_state, -1, 1);
    slap_MatMulAdd(Qx, slap_Transpose(ineq_state_jac), ineq_state2, 1, 1);
    // Qxx += ρ*∇hx'*mask*∇hx
    slap_MatMulAdd(ineq_state_jac2, mask_state, ineq_state_jac, 1, 0);
    slap_MatMulAdd(Qxx, slap_Transpose(ineq_state_jac), ineq_state_jac2, 1, 1);

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
                   1);  // P -= 2*K'Qux
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

  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob);
  // Replace P[N] since we used it for Quu_temp (need improving later)
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_MpcLti(Matrix* X, Matrix* U, tiny_ProblemData* prob,
                                tiny_Solver* solver, const tiny_LtiModel model,
                                const int verbose) {
  int N = prob->nhorizon;
  int n = prob->nstates;
  int m = prob->ninputs;
  for (int k = 0; k < N - 1; ++k) {
    tiny_DynamicsLti(&(X[k + 1]), X[k], U[k], model);
  }
  double G_temp_data[(n + m) * (n + m + 1)];
  Matrix Q_temp = slap_MatrixFromArray(n + m, n + m + 1, G_temp_data);

  double ineq_temp_data[prob->ncstr_states *
                        (prob->ncstr_states + prob->ncstr_states + 2)];
  Matrix ineq_temp = slap_MatrixFromArray(
      prob->ncstr_states, prob->ncstr_states + 2 * prob->nstates + 2,
      ineq_temp_data);

  Matrix ineq_input =
      slap_CreateSubMatrix(ineq_temp, 0, 0, prob->ncstr_inputs, 1);
  Matrix new_input_duals =
      slap_CreateSubMatrix(ineq_temp, 0, 1, prob->ncstr_inputs, 1);
  Matrix mask_input = slap_CreateSubMatrix(ineq_temp, 0, 2, prob->ncstr_inputs,
                                           prob->ncstr_inputs);

  Matrix ineq_state =
      slap_CreateSubMatrix(ineq_temp, 0, 0, prob->ncstr_states, 1);
  Matrix new_state_duals =
      slap_CreateSubMatrix(ineq_temp, 0, 1, prob->ncstr_states, 1);
  Matrix mask_state = slap_CreateSubMatrix(ineq_temp, 0, 2, prob->ncstr_states,
                                           prob->ncstr_states);

  Matrix eq_goal = slap_MatrixFromArray(prob->ncstr_goal, 1, ineq_temp_data);

  double cstr_violation = 0.0;
  for (int iter = 0; iter < solver->max_primal_iters; ++iter) {
    if (verbose > 1) printf("backward pass\n");
    tiny_ConstrainedBackwardPassLti(prob, *solver, model, X, U, &Q_temp,
                                    &ineq_temp);
    if (verbose > 1) printf("forward pass\n");
    tiny_ForwardPassLti(X, U, *prob, model);

    if (verbose > 0) {
      printf("iter     J           ΔJ         reg         ρ\n");
      printf("--------------------------------------------------\n");
      printf("%3d   %10.3e  %9.2e  %9.2e   %9.2e\n", iter, 0.0, 0.0,
             solver->regu, solver->penalty);
    }

    if (verbose > 1) printf("update duals and penalty\n");

    // For linear systems, only 1 iteration
    cstr_violation = 0.0;
    double norm_inf = 0.0;

    for (int k = 0; k < N - 1; ++k) {
      //========= Control constraints ==========
      tiny_IneqInputs(&ineq_input, *prob,
                      U[k]);  // ineq_input size = 2*NINPUTS
      tiny_ActiveIneqMask(&mask_input, prob->input_duals[k], ineq_input);
      slap_ScaleByConst(mask_input, solver->penalty);  // mask = ρ*mask
      // Constraint violation
      slap_ArgMax(ineq_input, &norm_inf);
      norm_inf = norm_inf > 0.0 ? norm_inf : 0.0;
      norm_inf = norm_inf * 2;
      // convio = max(convio,norm(hxv + abs.(hxv),Inf))
      cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
      // Update duals
      tiny_IneqInputsOffset(&ineq_input, *prob);  // g
      slap_MatrixCopy(new_input_duals, prob->input_duals[k]);
      slap_MatMulAdd(new_input_duals, mask_input, ineq_input, -1,
                     1);  //μ[k] - ρ*mask * g
      tiny_ClampIneqDuals(&(prob->input_duals)[k], new_input_duals);
    }

    for (int k = 0; k < N; ++k) {
      //========= State constraints ==========
      tiny_IneqStates(&ineq_state, *prob,
                      X[k]);  // ineq_input size = 2*NINPUTS
      tiny_ActiveIneqMask(&mask_state, prob->state_duals[k], ineq_state);
      slap_ScaleByConst(mask_state, solver->penalty);  // mask = ρ*mask
      // Constraint violation
      slap_ArgMax(ineq_state, &norm_inf);
      norm_inf = norm_inf > 0.0 ? norm_inf : 0.0;
      norm_inf = norm_inf * 2;
      // convio = max(convio,norm(hxv + abs.(hxv),Inf))
      cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
      // Update duals
      tiny_IneqStatesOffset(&ineq_state, *prob);  // g
      slap_MatrixCopy(new_state_duals, prob->state_duals[k]);
      slap_MatMulAdd(new_state_duals, mask_state, ineq_state, -1,
                     1);  //μ[k] - ρ*mask*g
      tiny_ClampIneqDuals(&(prob->state_duals)[k], new_state_duals);
    }

    //========= Goal constraints ==========
    slap_MatrixAddition(eq_goal, X[N - 1], prob->X_ref[N - 1], -1);
    norm_inf = slap_NormInf(eq_goal);
    cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
    // λ -= ρ*h
    slap_MatrixCopy(eq_goal, prob->X_ref[N - 1]);  // h = xg
    slap_MatrixAddition(prob->goal_dual, prob->goal_dual, eq_goal,
                        -solver->penalty);

    if (verbose > 0) printf("convio: %.6f \n\n", cstr_violation);
    if (cstr_violation < solver->cstr_tol) {
      if (verbose > 0) printf("SUCCESS!\n");
      solver->penalty = 1;  // reset penalty for next MPC
      return SLAP_NO_ERROR;
    }
    solver->penalty = solver->penalty * solver->penalty_mul;
  }
  return SLAP_NO_ERROR;
}