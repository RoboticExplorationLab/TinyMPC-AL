#include "mpc_ltv.h"

enum slap_ErrorCode tiny_ConstrainedBackwardPassLtv(
    tiny_ProblemData* prob, const tiny_Settings solver, const tiny_LtvModel model,
    const Matrix* X, const Matrix* U, Matrix* Q_temp, Matrix* c_temp) {
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

  Matrix cu = slap_CreateSubMatrix(*c_temp, 0, 0, 2 * m, 1);
  Matrix cu2 = slap_CreateSubMatrix(*c_temp, 0, 1, 2 * m, 1);
  Matrix cu_mask = slap_CreateSubMatrix(*c_temp, 0, 2, 2 * m, 2 * m);
  Matrix cu_jac =
      slap_CreateSubMatrix(*c_temp, 0, 2 * m + 2, 2 * m, m);
  Matrix cu_jac2 =
      slap_CreateSubMatrix(*c_temp, 0, 2 * m + 2 + m, 2 * m, m);

  Matrix cx = slap_CreateSubMatrix(*c_temp, 0, 0, 2 * n, 1);
  Matrix cx2 = slap_CreateSubMatrix(*c_temp, 0, 1, 2 * n, 1);
  Matrix cx_mask = slap_CreateSubMatrix(*c_temp, 0, 2, 2 * n, 2 * n);
  Matrix cx_jac =
      slap_CreateSubMatrix(*c_temp, 0, 2 * n + 2, 2 * n, n);
  Matrix cx_jac2 =
      slap_CreateSubMatrix(*c_temp, 0, 2 * n + 2 + n, 2 * n, n);

  // ========= Goal constraints ==========
  if (prob->ncstr_goal > 0) {
    slap_Copy(Qx, prob->X_ref[N - 1]);  // h = xg
    slap_SetIdentity(Qxx, 1);           // H constant
    slap_MatrixAddition(prob->p[N - 2], prob->YG, Qx,
                        -solver.penalty);  // (λ - ρ*h)
    slap_MatMulAdd(prob->p[N - 1], slap_Transpose(Qxx), prob->p[N - 2], 1,
                   1);  // p[N]  += H'*(λ - ρ*h)
    slap_MatMulAdd(prob->P[N - 1], slap_Transpose(Qxx), Qxx, solver.penalty,
                   1);  // P[N]  += ρ*H'H
  }

  //========= State constraints at end ==========
  if (prob->ncstr_states > 0) {
    tiny_EvalStateConstraint(&cx, *prob,
                    X[N - 1]);  // cx size = 2*NINPUTS
    tiny_ActiveIneqMask(&cx_mask, prob->YX[N - 1], cx);
    slap_ScaleByConst(cx_mask, solver.penalty);  // mask = ρ*mask
    // tiny_EvalStateConstraintJacobian(&cx_jac, *prob);
    // Qx  += G'*(μx[k] - ρ*mask * g)
    // tiny_EvalStateConstraintOffset(&cx, *prob);  // g
    slap_Copy(cx2, prob->YX[N - 1]);
    slap_MatMulAdd(cx2, cx_mask, prob->bcx, -1,
                   1);  //μx[k] - ρ*mask*g
    slap_MatMulAdd(prob->p[N - 1], slap_Transpose(prob->Acx),
                   cx2, 1, 1);
    // Qxx += G'*ρmask*G
    slap_MatMulAdd(cx_jac2, cx_mask, prob->Acx, 1, 0);
    slap_MatMulAdd(prob->P[N - 1], slap_Transpose(prob->Acx),
                   cx_jac2, 1, 1);
  }
  for (int k = N - 2; k >= 0; --k) {
    // Stage cost expansion
    tiny_ExpandStageCost(&Qxx, &Qx, &Quu, &Qu, *prob, k);
    // State Gradient: Qx = q + A'(P*f + p)
    slap_Copy(prob->p[k], prob->p[k + 1]);
    slap_MatMulAdd(prob->p[k], prob->P[k + 1], model.f[k], 1,
                   1);  // p[k] = P[k+1]*f + p[k+1]
    slap_MatMulAdd(Qx, slap_Transpose(model.A[k]), prob->p[k], 1, 1);

    // Control Gradient: Qu = r + B'(P*f + p)
    slap_MatMulAdd(Qu, slap_Transpose(model.B[k]), prob->p[k], 1, 1);

    // State Hessian: Qxx = Q + A'P*A
    slap_MatMulAdd(prob->P[k], prob->P[k + 1], model.A[k], 1,
                   0);  // P[k] = P[k+1]*A
    slap_MatMulAdd(Qxx, slap_Transpose(model.A[k]), prob->P[k], 1,
                   1);  // Qxx = Q + A'P*A

    // Control Hessian Quu = R + B'P*B
    slap_MatMulAdd(Qxu, prob->P[k + 1], model.B[k], 1, 0);  // Qxu = P * B
    slap_MatMulAdd(Quu, slap_Transpose(model.B[k]), Qxu, 1,
                   1);  // Quu = R + B'P*B
    slap_MatMulAdd(Quu, slap_Transpose(model.B[k]), model.B[k], solver.reg, 1);
    // Hessian Cross-Term
    slap_MatMulAdd(Qux, slap_Transpose(model.B[k]), prob->P[k], 1,
                   0);  // Qux = B'P*A

    //========= Control constraints ==========
    if (prob->ncstr_inputs > 0) {
      tiny_EvalInputConstraint(&cu, *prob, U[k]);  // cu size = 2*NINPUTS
      tiny_ActiveIneqMask(&cu_mask, prob->YU[k], cu);
      slap_ScaleByConst(cu_mask, solver.penalty);  // mask = ρ*mask
      // tiny_EvalInputConstraintJacobian(&cu_jac, *prob);
      // Qu  += G'*(μ[k] + (mask * g)
      // tiny_EvalInputConstraintOffset(&cu, *prob);  // g
      slap_Copy(cu2, prob->YU[k]);
      slap_MatMulAdd(cu2, cu_mask, prob->bcu, -1, 1);
      slap_MatMulAdd(Qu, slap_Transpose(prob->Acu), cu2, 1, 1);
      // Quu += ∇hu'*mask*∇hu
      slap_MatMulAdd(cu_jac2, cu_mask, prob->Acu, 1, 0);
      slap_MatMulAdd(Quu, slap_Transpose(prob->Acu), cu_jac2, 1,
                     1);
    }

    //========= State constraints ==========
    if (prob->ncstr_states > 0) {
      tiny_EvalStateConstraint(&cx, *prob, X[k]);  // cx size = 2*NINPUTS
      tiny_ActiveIneqMask(&cx_mask, prob->YX[k], cx);
      slap_ScaleByConst(cx_mask, solver.penalty);  // mask = ρ*mask
      // tiny_EvalStateConstraintJacobian(&cx_jac, *prob);
      // Qx  += G'*(μx[k] - ρ*mask * g)
      // tiny_EvalStateConstraintOffset(&cx, *prob);  // g
      slap_Copy(cx2, prob->YX[k]);
      slap_MatMulAdd(cx2, cx_mask, prob->bcx, -1, 1);
      slap_MatMulAdd(Qx, slap_Transpose(prob->Acx), cx2, 1, 1);
      // Qxx += ρ*∇hx'*mask*∇hx
      slap_MatMulAdd(cx_jac2, cx_mask, prob->Acx, 1, 0);
      slap_MatMulAdd(Qxx, slap_Transpose(prob->Acx), cx_jac2, 1,
                     1);
    }

    // Calculate Gains
    slap_Copy(Quu_temp, Quu);
    slap_Cholesky(Quu_temp);
    slap_Copy(prob->K[k], Qux);
    slap_Copy(prob->d[k], Qu);
    slap_CholeskySolve(Quu_temp, prob->d[k]);  // d = Quu\Qu
    slap_CholeskySolve(Quu_temp, prob->K[k]);  // K = Quu\Qux

    // Cost-to-Go Hessian: P = Qxx + K'Quu*K - K'Qux - Qux'K
    slap_Copy(prob->P[k], Qxx);                                  // P = Qxx
    slap_MatMulAdd(Qxu, slap_Transpose(prob->K[k]), Quu, 1, 0);  // Qxu = K'Quu
    slap_MatMulAdd(prob->P[k], Qxu, prob->K[k], 1, 1);           // P += K'Quu*K
    slap_MatMulAdd(prob->P[k], slap_Transpose(prob->K[k]), Qux, -2,
                   1);  // P -= 2*K'Qux
    // slap_MatMulAdd(prob->P[k], slap_Transpose(Qux), prob->K[k], -1,
    //                1);  // P -= Qux'K

    // Cost-to-Go Gradient: p = Qx + K'Quu*d - K'Qu - Qux'd
    slap_Copy(prob->p[k], Qx);                          // p = Qx
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

enum slap_ErrorCode tiny_MpcLtv(Matrix* X, Matrix* U, tiny_ProblemData* prob,
                                tiny_Settings* solver, const tiny_LtvModel model,
                                const int verbose, sfloat* temp_data) {
  int N = prob->nhorizon;
  int n = prob->nstates;
  int m = prob->ninputs;

  for (int k = 0; k < N - 1; ++k) {
    tiny_DynamicsLtv(&(X[k + 1]), X[k], U[k], model, k);
  }

  // sfloat temp_data[(n + m) * (n + m + 1)];
  Matrix Q_temp = slap_MatrixFromArray(n + m, n + m + 1, temp_data);

  // sfloat ineq_temp_data[2*n *
  // (2*n + 2*n + 2)];
  Matrix c_temp = slap_MatrixFromArray(2 * n, 2 * n + 2 * n + 2,
                                          &temp_data[(n + m) * (n + m + 1)]);

  Matrix cu = slap_CreateSubMatrix(c_temp, 0, 0, 2 * m, 1);
  Matrix YU_hat = slap_CreateSubMatrix(c_temp, 0, 1, 2 * m, 1);
  Matrix cu_mask = slap_CreateSubMatrix(c_temp, 0, 2, 2 * m, 2 * m);

  Matrix cx = slap_CreateSubMatrix(c_temp, 0, 0, 2 * n, 1);
  Matrix YX_hat = slap_CreateSubMatrix(c_temp, 0, 1, 2 * n, 1);
  Matrix cx_mask = slap_CreateSubMatrix(c_temp, 0, 2, 2 * n, 2 * n);

  Matrix cg =
      slap_MatrixFromArray(n, 1, &temp_data[2 * n * (2 * n + 2 * n + 2)]);

  sfloat cstr_violation = 0.0;
  for (int iter = 0; iter < solver->max_outer_iters; ++iter) {
    if (verbose > 1) printf("backward pass\n");
    tiny_ConstrainedBackwardPassLtv(prob, *solver, model, X, U, &Q_temp,
                                    &c_temp);
    if (verbose > 1) printf("forward pass\n");
    tiny_ForwardPassLtv(X, U, *prob, model);

    if (verbose > 1) printf("update duals and penalty\n");

    // For linear systems, only 1 iteration
    cstr_violation = 0.0;
    sfloat norm_inf = 0.0;

    if (prob->ncstr_inputs > 0) {
      for (int k = 0; k < N - 1; ++k) {
        //========= Control constraints ==========
        tiny_EvalInputConstraint(&cu, *prob,
                        U[k]);  // cu size = 2*NINPUTS
        tiny_ActiveIneqMask(&cu_mask, prob->YU[k], cu);
        slap_ScaleByConst(cu_mask, solver->penalty);  // mask = ρ*mask
        // Constraint violation
        slap_ArgMax(cu, &norm_inf);
        norm_inf = norm_inf > 0.0 ? norm_inf : 0.0;
        norm_inf = norm_inf * 2;
        // convio = max(convio,norm(hxv + abs.(hxv),Inf))
        cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
        // Update duals
        // tiny_EvalInputConstraintOffset(&cu, *prob);  // g
        slap_Copy(YU_hat, prob->YU[k]);
        slap_MatMulAdd(YU_hat, cu_mask, prob->bcu, -1,
                       1);  //μ[k] - ρ*mask * g
        tiny_ProjectOrthantDuals(&(prob->YU)[k], YU_hat);
      }
    }

    if (prob->ncstr_states > 0) {
      for (int k = 0; k < N; ++k) {
        //========= State constraints ==========
        tiny_EvalStateConstraint(&cx, *prob,
                        X[k]);  // cu size = 2*NINPUTS
        tiny_ActiveIneqMask(&cx_mask, prob->YX[k], cx);
        slap_ScaleByConst(cx_mask, solver->penalty);  // mask = ρ*mask
        // Constraint violation
        slap_ArgMax(cx, &norm_inf);
        norm_inf = norm_inf > 0.0 ? norm_inf : 0.0;
        norm_inf = norm_inf * 2;
        // convio = max(convio,norm(hxv + abs.(hxv),Inf))
        cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
        // Update duals
        // tiny_EvalStateConstraintOffset(&cx, *prob);  // g
        slap_Copy(YX_hat, prob->YX[k]);
        slap_MatMulAdd(YX_hat, cx_mask, prob->bcx, -1,
                       1);  //μ[k] - ρ*mask*g
        tiny_ProjectOrthantDuals(&(prob->YX)[k], YX_hat);
      }
    }

    if (prob->ncstr_goal > 0) {
      // ========= Goal constraints ==========
      slap_MatrixAddition(cg, X[N - 1], prob->X_ref[N - 1], -1);
      norm_inf = slap_NormInf(cg);
      cstr_violation = cstr_violation < norm_inf ? norm_inf : cstr_violation;
      // λ -= ρ*h
      slap_Copy(cg, prob->X_ref[N - 1]);  // h = xg
      slap_MatrixAddition(prob->YG, prob->YG, cg,
                          -solver->penalty);
    }
    if (verbose > 0) {
      if (iter == 0) {
        printf("iter     J           ΔJ         convio        reg         ρ\n");
        printf("-----------------------------------------------------------\n");
      }
      printf("%3d   %10.3e  %9.2e  %9.2e %10.1e  %8.1e\n", iter, 0.0, 0.0,
             cstr_violation, solver->reg, solver->penalty);
    }
    if (cstr_violation < solver->cstr_tol) {
      if (verbose > 0) printf("SUCCESS!\n");
      solver->penalty = 1;  // reset penalty for next MPC
      return SLAP_NO_ERROR;
    }
    solver->penalty = solver->penalty * solver->penalty_mul;
  }
  solver->penalty = 1;  // reset penalty for next MPC
  return SLAP_NO_ERROR;
}