#include "nmpc.h"

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
                          const tiny_ProblemData prob, const Matrix x,
                          const Matrix u, const int k) {
  double dx_data[prob.nstates];
  Matrix dx = slap_MatrixFromArray(prob.nstates, 1, dx_data);
  slap_MatrixAddition(dx, x, prob.X_ref[k], -1);
  slap_Copy(*hes_el_xx, prob.Q);
  slap_MatMulAdd(*grad_el_x, prob.Q, dx, 1, 0);
  Matrix du = slap_MatrixFromArray(prob.ninputs, 1, dx_data);
  slap_MatrixAddition(du, u, prob.U_ref[k], -1);
  slap_Copy(*hes_el_uu, prob.R);
  slap_MatMulAdd(*grad_el_u, prob.R, du, 1, 0);
}

void tiny_ExpandTerminalCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                             const tiny_ProblemData prob, const Matrix x) {
  double dx_data[prob.nstates];
  Matrix dx = slap_MatrixFromArray(prob.nstates, 1, dx_data);
  slap_MatrixAddition(dx, x, prob.X_ref[prob.nhorizon - 1], -1);
  slap_Copy(*hes_el_xx, prob.Qf);
  slap_MatMulAdd(*grad_el_x, prob.Qf, dx, 1, 0);
}

// Riccati recursion for LTI without constraints
enum slap_ErrorCode tiny_BackwardPassLti(tiny_ProblemData* prob,
                                         const tiny_LtiModel model,
                                         const tiny_Solver solver,
                                         const Matrix* X, const Matrix* U,
                                         Matrix G_temp) {
  // Copy terminal cost-to-go
  int N = prob->nhorizon;
  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob,
                          X[N - 1]);
  int n = prob->nstates;
  int m = prob->ninputs;

  Matrix Gxx = slap_CreateSubMatrix(G_temp, 0, 0, n, n);
  Matrix Gxu = slap_CreateSubMatrix(G_temp, 0, n, n, m);
  Matrix Gux = slap_CreateSubMatrix(G_temp, n, 0, m, n);
  Matrix Guu = slap_CreateSubMatrix(G_temp, n, n, m, m);
  Matrix Gx = slap_CreateSubMatrix(G_temp, 0, n + m, n, 1);
  Matrix Gu = slap_CreateSubMatrix(G_temp, n, n + m, m, 1);

  // Hijack the first part of P[N] to use for Cholesky decomposition
  // NOTE: Assumes m <= n
  Matrix Quu_temp = slap_Reshape(prob->P[N - 1], m, m);
  for (int k = N - 2; k >= 0; --k) {
    // Stage cost expansion
    tiny_ExpandStageCost(&Gxx, &Gx, &Guu, &Gu, *prob, X[k], U[k], k);
    // State Gradient: Gx = q + A'(P*f + p)
    slap_Copy(prob->p[k], prob->p[k + 1]);
    // slap_MatMulAdd(prob->p[k], prob->P[k + 1], model.f, 1,
    //  1);  // p[k] = P[k+1]*f + p[k+1]
    slap_MatMulAdd(Gx, slap_Transpose(model.A), prob->p[k], 1, 1);

    // Control Gradient: Gu = r + B'(P*f + p)
    slap_MatMulAdd(Gu, slap_Transpose(model.B), prob->p[k], 1, 1);

    // State Hessian: Gxx = Q + A'P*A
    slap_MatMulAdd(prob->P[k], prob->P[k + 1], model.A, 1,
                   0);  // P[k] = P[k+1]*A
    slap_MatMulAdd(Gxx, slap_Transpose(model.A), prob->P[k], 1,
                   1);  // Gxx = Q + A'P*A

    // Control Hessian Guu = R + B'P*B
    slap_MatMulAdd(Gxu, prob->P[k + 1], model.B, 1, 0);       // Gxu = P * B
    slap_MatMulAdd(Guu, slap_Transpose(model.B), Gxu, 1, 1);  // Guu = R + B'P*B
    slap_MatMulAdd(Guu, slap_Transpose(model.B), model.B, solver.reg, 1);
    // Hessian Cross-Term
    slap_MatMulAdd(Gux, slap_Transpose(model.B), prob->P[k], 1,
                   0);  // Gux = B'P*A

    // Calculate Gains
    slap_Copy(Quu_temp, Guu);
    slap_Cholesky(Quu_temp);
    slap_Copy(prob->K[k], Gux);
    slap_Copy(prob->d[k], Gu);
    slap_CholeskySolve(Quu_temp, prob->d[k]);  // d = Guu\Gu
    slap_CholeskySolve(Quu_temp, prob->K[k]);  // K = Guu\Gux

    // Cost-to-Go Hessian: P = Gxx + K'Guu*K - K'Gux - Gux'K
    slap_Copy(prob->P[k], Gxx);                            // P = Gxx
    slap_MatMulAdd(Gxu, slap_Transpose(prob->K[k]), Guu, 1, 0);  // Gxu = K'Guu
    slap_MatMulAdd(prob->P[k], Gxu, prob->K[k], 1, 1);           // P += K'Guu*K
    slap_MatMulAdd(prob->P[k], slap_Transpose(prob->K[k]), Gux, -2,
                   1);  // P -= K'Gux
    // slap_MatMulAdd(prob->P[k], slap_Transpose(Gux), prob->K[k], -1,
    //                1);  // P -= Gux'K

    // Cost-to-Go Gradient: p = Gx + K'Guu*d - K'Gu - Gux'd
    slap_Copy(prob->p[k], Gx);                    // p = Gx
    slap_MatMulAdd(prob->p[k], Gxu, prob->d[k], 1, 1);  // p += K'Guu*d
    slap_MatMulAdd(prob->p[k], slap_Transpose(prob->K[k]), Gu, -1,
                   1);  // p -= K'Gu
    slap_MatMulAdd(prob->p[k], slap_Transpose(Gux), prob->d[k], -1,
                   1);  // p -= Gux'd
  }
  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob,
                          X[N - 1]);  // Replace P[N] since we used it for
                                      // Quu_temp (need improving later)
  return SLAP_NO_ERROR;
}

// Riccati recursion for LTI with constraints
enum slap_ErrorCode tiny_ConstrainedBackwardPassLti(
    tiny_ProblemData* prob, const tiny_LtiModel model, const tiny_Solver solver,
    const Matrix* X, const Matrix* U, Matrix* G_temp, Matrix* ineq_temp) {
  // Copy terminal cost-to-go
  int N = prob->nhorizon;
  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob,
                          X[N - 1]);
  int n = prob->nstates;
  int m = prob->ninputs;

  Matrix Gxx = slap_CreateSubMatrix(*G_temp, 0, 0, n, n);
  Matrix Gxu = slap_CreateSubMatrix(*G_temp, 0, n, n, m);
  Matrix Gux = slap_CreateSubMatrix(*G_temp, n, 0, m, n);
  Matrix Guu = slap_CreateSubMatrix(*G_temp, n, n, m, m);
  Matrix Gx = slap_CreateSubMatrix(*G_temp, 0, n + m, n, 1);
  Matrix Gu = slap_CreateSubMatrix(*G_temp, n, n + m, m, 1);

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
  slap_MatrixAddition(Gx, X[N - 1], prob->X_ref[N - 1], -1);  // hxv
  slap_SetIdentity(Gxx, 1);                                   // grad_hx
  slap_MatrixAddition(prob->p[N - 2], prob->goal_dual, Gx,
                      solver.penalty);  // p[N]  += ∇hx'*(λ + ρ*hxv)
  slap_MatMulAdd(prob->p[N - 1], slap_Transpose(Gxx), prob->p[N - 2], 1,
                 1);  // !!constant jacobian
  slap_MatMulAdd(prob->P[N - 1], slap_Transpose(Gxx), Gxx, solver.penalty,
                 1);  // P[N]  += ρ*∇hx'∇hx

  //========= State constraints at end ==========
  tiny_IneqStates(&ineq_state, *prob, X[N - 1]);  // ineq_state size = 2*NINPUTS
  tiny_ActiveIneqMask(&mask_state, prob->state_duals[N - 1], ineq_state);
  slap_ScaleByConst(mask_state, solver.penalty);  // mask = ρ*mask
  tiny_IneqStatesJacobian(&ineq_state_jac, *prob);
  // Gx  += ∇hx'*(μx[k] + ρ*(mask * hxv))
  slap_Copy(ineq_state2, prob->state_duals[N - 1]);
  slap_MatMulAdd(ineq_state2, mask_state, ineq_state, 1, 1);
  slap_MatMulAdd(prob->p[N - 1], slap_Transpose(ineq_state_jac), ineq_state2, 1,
                 1);
  // Gxx += ρ*∇hx'*mask*∇hx
  slap_MatMulAdd(ineq_state_jac2, mask_state, ineq_state_jac, 1, 0);
  slap_MatMulAdd(prob->P[N - 1], slap_Transpose(ineq_state_jac),
                 ineq_state_jac2, 1, 1);

  for (int k = N - 2; k >= 0; --k) {
    // Stage cost expansion
    tiny_ExpandStageCost(&Gxx, &Gx, &Guu, &Gu, *prob, X[k], U[k], k);
    // State Gradient: Gx = q + A'(P*f + p)
    slap_Copy(prob->p[k], prob->p[k + 1]);
    // slap_MatMulAdd(prob->p[k], prob->P[k + 1], model.f, 1,
    //                1);  // p[k] = P[k+1]*f + p[k+1]
    slap_MatMulAdd(Gx, slap_Transpose(model.A), prob->p[k], 1, 1);

    // Control Gradient: Gu = r + B'(P*f ineq_input_jac+ p)
    slap_MatMulAdd(Gu, slap_Transpose(model.B), prob->p[k], 1, 1);

    // State Hessian: Gxx = Q + A'P*A
    slap_MatMulAdd(prob->P[k], prob->P[k + 1], model.A, 1,
                   0);  // P[k] = P[k+1]*A
    slap_MatMulAdd(Gxx, slap_Transpose(model.A), prob->P[k], 1,
                   1);  // Gxx = Q + A'P*A

    // Control Hessian Guu = R + B'P*B
    slap_MatMulAdd(Gxu, prob->P[k + 1], model.B, 1, 0);       // Gxu = P * B
    slap_MatMulAdd(Guu, slap_Transpose(model.B), Gxu, 1, 1);  // Guu = R + B'P*B
    slap_MatMulAdd(Guu, slap_Transpose(model.B), model.B, solver.reg, 1);
    // Hessian Cross-Term
    slap_MatMulAdd(Gux, slap_Transpose(model.B), prob->P[k], 1,
                   0);  // Gux = B'P*A

    //========= Control constraints ==========
    tiny_IneqInputs(&ineq_input, *prob, U[k]);  // ineq_input size = 2*NINPUTS
    tiny_ActiveIneqMask(&mask_input, prob->input_duals[k], ineq_input);
    slap_ScaleByConst(mask_input, solver.penalty);  // mask = ρ*mask
    tiny_IneqInputsJacobian(&ineq_input_jac, *prob);
    // Gu  += ∇hu'*(μ[k] + (mask * huv))
    slap_Copy(ineq_input2, prob->input_duals[k]);
    slap_MatMulAdd(ineq_input2, mask_input, ineq_input, 1, 1);
    slap_MatMulAdd(Gu, slap_Transpose(ineq_input_jac), ineq_input2, 1, 1);
    // Guu += ∇hu'*mask*∇hu
    slap_MatMulAdd(ineq_input_jac2, mask_input, ineq_input_jac, 1, 0);
    slap_MatMulAdd(Guu, slap_Transpose(ineq_input_jac), ineq_input_jac2, 1, 1);

    //========= State constraints ==========
    tiny_IneqStates(&ineq_state, *prob, X[k]);  // ineq_state size = 2*NINPUTS
    tiny_ActiveIneqMask(&mask_state, prob->state_duals[k], ineq_state);
    slap_ScaleByConst(mask_state, solver.penalty);  // mask = ρ*mask
    tiny_IneqStatesJacobian(&ineq_state_jac, *prob);
    // Gx  += ∇hx'*(μx[k] + ρ*(mask * hxv))
    slap_Copy(ineq_state2, prob->state_duals[k]);
    slap_MatMulAdd(ineq_state2, mask_state, ineq_state, 1, 1);
    slap_MatMulAdd(Gx, slap_Transpose(ineq_state_jac), ineq_state2, 1, 1);
    // Gxx += ρ*∇hx'*mask*∇hx
    slap_MatMulAdd(ineq_state_jac2, mask_state, ineq_state_jac, 1, 0);
    slap_MatMulAdd(Gxx, slap_Transpose(ineq_state_jac), ineq_state_jac2, 1, 1);

    // Calculate Gains
    slap_Copy(Quu_temp, Guu);
    slap_Cholesky(Quu_temp);
    slap_Copy(prob->K[k], Gux);
    slap_Copy(prob->d[k], Gu);
    slap_CholeskySolve(Quu_temp, prob->d[k]);  // d = Guu\Gu
    slap_CholeskySolve(Quu_temp, prob->K[k]);  // K = Guu\Gux

    // Cost-to-Go Hessian: P = Gxx + K'Guu*K - K'Gux - Gux'K
    slap_Copy(prob->P[k], Gxx);                            // P = Gxx
    slap_MatMulAdd(Gxu, slap_Transpose(prob->K[k]), Guu, 1, 0);  // Gxu = K'Guu
    slap_MatMulAdd(prob->P[k], Gxu, prob->K[k], 1, 1);           // P += K'Guu*K
    slap_MatMulAdd(prob->P[k], slap_Transpose(prob->K[k]), Gux, -2,
                   1);  // P -= 2*K'Gux
    // slap_MatMulAdd(prob->P[k], slap_Transpose(Gux), prob->K[k], -1,
    //                1);  // P -= Gux'K

    // Cost-to-Go Gradient: p = Gx + K'Guu*d - K'Gu - Gux'd
    slap_Copy(prob->p[k], Gx);                    // p = Gx
    slap_MatMulAdd(prob->p[k], Gxu, prob->d[k], 1, 1);  // p += K'Guu*d
    slap_MatMulAdd(prob->p[k], slap_Transpose(prob->K[k]), Gu, -1,
                   1);  // p -= K'Gu
    slap_MatMulAdd(prob->p[k], slap_Transpose(Gux), prob->d[k], -1,
                   1);  // p -= Gux'd
  }

  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob,
                          X[N - 1]);  // Replace P[N] since we used it for
                                      // Quu_temp (need improving later)
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_AugmentedLagrangianLqr(Matrix* X, Matrix* U,
                                                tiny_ProblemData* prob,
                                                tiny_Solver* solver,
                                                const tiny_LtiModel model,
                                                const int verbose) {
  int N = prob->nhorizon;
  int n = prob->nstates;
  int m = prob->ninputs;
  for (int k = 0; k < N - 1; ++k) {
    tiny_DynamicsLti(&(X[k + 1]), X[k], U[k], model);
  }
  double G_temp_data[(n + m) * (n + m + 1)];
  Matrix G_temp = slap_MatrixFromArray(n + m, n + m + 1, G_temp_data);

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

  double norm_d_max = 0.0;
  double cstr_violation = 0.0;
  for (int iter = 0; iter < solver->max_primal_iters; ++iter) {
    tiny_ConstrainedBackwardPassLti(prob, model, *solver, X, U, &G_temp,
                                    &ineq_temp);
    tiny_ForwardPassLti(X, U, *prob, model);

    norm_d_max = tiny_RiccatiConvergence(*prob);

    if (verbose == 1) {
      printf("outer loop\n");
      printf(
          "iter     J           ΔJ        |d|         α        reg         "
          "ρ\n");
      printf(
          "--------------------------------------------------------------------"
          "-\n");
      printf("%3d   %10.3e  %9.2e  %9.2e  %6.4f   %9.2e   %9.2e\n", iter, 0.0,
             0.0, norm_d_max, 1.0, solver->reg, solver->penalty);
    }

    printf("update duals and penalty\n");

    // For linear systems, only 1 iteration, shouldn't need condition here
    if (0 * norm_d_max < solver->riccati_tol) {
      cstr_violation = 0.0;
      double norm_inf = 0.0;

      for (int k = 0; k < N - 1; ++k) {
        //========= Control constraints ==========
        tiny_IneqInputs(&ineq_input, *prob,
                        U[k]);  // ineq_input size = 2*NINPUTS
        tiny_ActiveIneqMask(&mask_input, prob->input_duals[k], ineq_input);
        slap_ScaleByConst(mask_input, solver->penalty);  // mask = ρ*mask
        // Update duals
        slap_Copy(new_input_duals, prob->input_duals[k]);
        slap_MatMulAdd(new_input_duals, mask_input, ineq_input, 1,
                       1);  //μ[k] + ρ*mask*huv
        slap_ArgMax(ineq_input, &norm_inf);
        norm_inf = norm_inf > 0.0 ? norm_inf : 0.0;
        norm_inf = norm_inf * 2;
        tiny_ClampIneqDuals(&(prob->input_duals)[k], new_input_duals);
        // convio = max(convio,norm(hxv + abs.(hxv),Inf))
        cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
      }

      for (int k = 0; k < N; ++k) {
        //========= State constraints ==========
        tiny_IneqStates(&ineq_state, *prob,
                        X[k]);  // ineq_input size = 2*NINPUTS
        tiny_ActiveIneqMask(&mask_state, prob->state_duals[k], ineq_state);
        slap_ScaleByConst(mask_state, solver->penalty);  // mask = ρ*mask
        // Update duals
        slap_Copy(new_state_duals, prob->state_duals[k]);
        slap_MatMulAdd(new_state_duals, mask_state, ineq_state, 1,
                       1);  //μ[k] + ρ*mask*huv
        slap_ArgMax(ineq_state, &norm_inf);
        norm_inf = norm_inf > 0.0 ? norm_inf : 0.0;
        norm_inf = norm_inf * 2;
        tiny_ClampIneqDuals(&(prob->state_duals)[k], new_state_duals);
        // convio = max(convio,norm(hxv + abs.(hxv),Inf))
        cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
      }

      //========= Goal constraints ==========
      slap_MatrixAddition(eq_goal, X[N - 1], prob->X_ref[N - 1], -1);
      norm_inf = slap_NormInf(eq_goal);
      cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
      // λ .+= ρ*hxv
      slap_MatrixAddition(prob->goal_dual, prob->goal_dual, eq_goal,
                          solver->penalty);

      printf("convio: %.6f \n", cstr_violation);
      if (cstr_violation < solver->cstr_tol) {
        printf("\nSUCCESS!\n");
        return SLAP_NO_ERROR;
      }
      solver->penalty = solver->penalty * solver->penalty_mul;
    }
  }
  return SLAP_NO_ERROR;
}

// Roll out the closed-loop dynamics with K and d from backward pass,
// calculate new X, U in place
enum slap_ErrorCode tiny_ForwardPassLti(Matrix* X, Matrix* U,
                                        const tiny_ProblemData prob,
                                        const tiny_LtiModel model) {
  int N = prob.nhorizon;
  double u_data[prob.ninputs];
  double x_data[prob.nstates];
  Matrix Un = slap_MatrixFromArray(prob.ninputs, 1, u_data);
  Matrix Xn = slap_MatrixFromArray(prob.nstates, 1, x_data);
  slap_Copy(Xn, prob.x0);
  for (int k = 0; k < N - 1; ++k) {
    // delta_x and delta_u over previous X, U
    // Control input: u = uf - d - K*(x - xf)
    slap_MatrixAddition(Un, U[k], prob.d[k], -1);  // u[k] = uf - d[k]
    slap_MatMulAdd(Un, prob.K[k], X[k], -1, 1);    // u[k] -= K[k] * x[k]
    slap_MatMulAdd(Un, prob.K[k], Xn, 1, 1);       // u[k] += K[k] * xf
    slap_Copy(U[k], Un);
    // Next state: x = A*x + B*u + f
    slap_Copy(Xn, X[k + 1]);
    tiny_DynamicsLti(&X[k + 1], X[k], U[k], model);
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

void tiny_DynamicsLti(Matrix* xn, const Matrix x, const Matrix u,
                      const tiny_LtiModel model) {
  slap_MatMulAdd(*xn, model.A, x, 1, 0);      // x[k+1] = A * x[k]
  slap_MatMulAdd(*xn, model.B, u, 1, 1);      // x[k+1] += B * u[k]
  slap_MatrixAddition(*xn, *xn, model.f, 1);  // x[k+1] += f
}