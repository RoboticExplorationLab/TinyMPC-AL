#include "al_lqr.h"

enum tiny_ErrorCode tiny_ConstrainedForwardPass(tiny_Workspace* work) {
  // tiny_RollOutClosedLoop(work);
  tiny_RollOutClosedLoopCost(work);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_ConstrainedBackwardPass(tiny_Workspace* work) {
  tiny_Model* model = work->data->model;
  int N = model[0].nhorizon;
  // int n = model[0].nstates;
  int m = model[0].ninputs;

  tiny_ExpandTerminalCost(work);
  // Hijack the first part of P[N] to use for Cholesky decomposition
  // NOTE: Assumes m <= n
  Matrix Quu_temp = slap_Reshape(work->soln->P[N - 1], m, m);

  //========= Goal constraints ==========
  if (work->stgs->en_cstr_goal) {
    slap_Copy(work->Qx, work->data->X_ref[N - 1]);  // h = xg
    slap_SetIdentity(work->Qxx, 1);           // H constant
    slap_MatrixAddition(work->soln->p[N - 2], work->soln->YG, work->Qx, -work->penalty);  // (λ - ρ*h)
    slap_MatMulAdd(work->soln->p[N - 1], slap_Transpose(work->Qxx), work->soln->p[N - 2], 1, 1);  // p[N]  += H'*(λ - ρ*h)
    slap_MatMulAdd(work->soln->P[N - 1], slap_Transpose(work->Qxx), work->Qxx, work->penalty, 1);  // P[N]  += ρ*H'H
  }

  //========= State constraints at end ==========
  if (work->stgs->en_cstr_states) {
    tiny_EvalStateConstraint(work, N - 1);  // cx size = 2*NINPUTS
    tiny_ActiveIneqMask(&(work->cx_mask), work->soln->YX[N - 1], work->cx);
    slap_ScaleByConst(work->cx_mask, work->penalty);  // mask = ρ*mask
    // tiny_EvalStateConstraintJacobian(&cx_jac, *prob);
    // Qx  += G'*(μx[k] - ρ*mask * g)
    // tiny_EvalStateConstraintOffset(&work->cx, *prob);  // g
    slap_Copy(work->cx2, work->soln->YX[N - 1]);
    slap_MatMulAdd(work->cx2, work->cx_mask, work->data->bcx, -1, 1);  //μx[k] - ρ*mask*g
    slap_MatMulAdd(work->soln->p[N - 1], slap_Transpose(work->data->Acx), work->cx2, 1, 1);
    // Qxx += G'*ρmask*G
    slap_MatMulAdd(work->cx_jac2, work->cx_mask, work->data->Acx, 1, 0);
    slap_MatMulAdd(work->soln->P[N - 1], slap_Transpose(work->data->Acx), work->cx_jac2, 1, 1);
  }

  if (model[0].ltv && model[0].affine) {
    for (int k = N - 2; k >= 0; --k) {
      // Stage cost expansion: Qxx = Q; Qx = q; Quu = R; Qu = r;
      tiny_ExpandStageCost(work, k);
      // State Gradient: Qx = q + A'(P*f + p)
      slap_Copy(work->soln->p[k], work->soln->p[k + 1]);
      slap_MatMulAdd(work->soln->p[k], work->soln->P[k + 1], model[0].f[k], 1,
                    1);  // p[k] = P[k+1]*f + p[k+1]
      slap_MatMulAdd(work->Qx, slap_Transpose(model[0].A[k]), work->soln->p[k], 1, 1);

      // Control Gradient: Qu = r + B'(P*f + p)
      slap_MatMulAdd(work->Qu, slap_Transpose(model[0].B[k]), work->soln->p[k], 1, 1);

      // State Hessian: Qxx = Q + A'P*A
      slap_MatMulAdd(work->soln->P[k], work->soln->P[k + 1], model[0].A[k], 1,
                    0);  // P[k] = P[k+1]*A
      slap_MatMulAdd(work->Qxx, slap_Transpose(model[0].A[k]), work->soln->P[k], 1,
                    1);  // Qxx = Q + A'P*A

      // Control Hessian Quu = R + B'P*B
      slap_MatMulAdd(work->Qxu, work->soln->P[k + 1], model[0].B[k], 1, 0);       // Qxu = P * B
      slap_MatMulAdd(work->Quu, slap_Transpose(model[0].B[k]), work->Qxu, 1, 1);  // Quu = R + B'P*B
      slap_MatMulAdd(work->Quu, slap_Transpose(model[0].B[k]), model[0].B[k], work->reg, 1);
      // Hessian Cross-Term
      slap_MatMulAdd(work->Qux, slap_Transpose(model[0].B[k]), work->soln->P[k], 1,
                    0);  // Qux = B'P*A

      //========= Control constraints ==========
      if (work->stgs->en_cstr_inputs) {
        tiny_EvalInputConstraint(work, k);  // work->cu size = 2*NINPUTS
        tiny_ActiveIneqMask(&(work->cu_mask), work->soln->YU[k], work->cu);
        slap_ScaleByConst(work->cu_mask, work->penalty);  // mask = ρ*mask
        // tiny_EvalInputConstraintJacobian(&cu_jac, *prob);
        // Qu  += G'*(μ[k] + (mask * g)
        // tiny_EvalInputConstraintOffset(&work->cu, *prob);  // g
        slap_Copy(work->cu2, work->soln->YU[k]);
        slap_MatMulAdd(work->cu2, work->cu_mask, work->data->bcu, -1, 1);
        slap_MatMulAdd(work->Qu, slap_Transpose(work->data->Acu), work->cu2, 1, 1);
        // Quu += ∇hu'*mask*∇hu
        slap_MatMulAdd(work->cu_jac2, work->cu_mask, work->data->Acu, 1, 0);
        slap_MatMulAdd(work->Quu, slap_Transpose(work->data->Acu), work->cu_jac2, 1, 1);
      }

      //========= State constraints ==========
      if (work->stgs->en_cstr_states) {
        tiny_EvalStateConstraint(work, k);  // cx size = 2*NINPUTS
        tiny_ActiveIneqMask(&(work->cx_mask), work->soln->YX[k], work->cx);
        slap_ScaleByConst(work->cx_mask, work->penalty);  // mask = ρ*mask
        // tiny_EvalStateConstraintJacobian(&cx_jac, *prob);
        // Qx  += G'*(μx[k] - ρ*mask * g)
        // tiny_EvalStateConstraintOffset(&cx, *prob);  // g
        slap_Copy(work->cx2, work->soln->YX[k]);
        slap_MatMulAdd(work->cx2, work->cx_mask, work->data->bcx, -1, 1);
        slap_MatMulAdd(work->Qx, slap_Transpose(work->data->Acx), work->cx2, 1, 1);
        // Qxx += ρ*∇hx'*mask*∇hx
        slap_MatMulAdd(work->cx_jac2, work->cx_mask, work->data->Acx, 1, 0);
        slap_MatMulAdd(work->Qxx, slap_Transpose(work->data->Acx), work->cx_jac2, 1, 1);
      }

      // Calculate Gains
      slap_Copy(Quu_temp, work->Quu);
      slap_Cholesky(Quu_temp);
      slap_Copy(work->soln->K[k], work->Qux);
      slap_Copy(work->soln->d[k], work->Qu);
      slap_CholeskySolve(Quu_temp, work->soln->d[k]);  // d = Quu\Qu
      slap_CholeskySolve(Quu_temp, work->soln->K[k]);  // K = Quu\Qux

      // Cost-to-Go Hessian: P = Qxx + K'Quu*K - K'Qux - Qux'K
      slap_Copy(work->soln->P[k], work->Qxx);                                  // P = Qxx
      slap_MatMulAdd(work->Qxu, slap_Transpose(work->soln->K[k]), work->Quu, 1, 0);  // Qxu = K'Quu
      slap_MatMulAdd(work->soln->P[k], work->Qxu, work->soln->K[k], 1, 1);           // P += K'Quu*K
      slap_MatMulAdd(work->soln->P[k], slap_Transpose(work->soln->K[k]), work->Qux, -2,
                    1);  // P -= K'Qux
      // slap_MatMulAdd(work->soln->P[k], slap_Transpose(Qux), work->soln->K[k], -1,
      //                1);  // P -= Qux'K

      // Cost-to-Go Gradient: p = Qx + K'Quu*d - K'Qu - Qux'd
      slap_Copy(work->soln->p[k], work->Qx);                          // p = Qx
      slap_MatMulAdd(work->soln->p[k], work->Qxu, work->soln->d[k], 1, 1);  // p += K'Quu*d
      slap_MatMulAdd(work->soln->p[k], slap_Transpose(work->soln->K[k]), work->Qu, -1,
                    1);  // p -= K'Qu
      slap_MatMulAdd(work->soln->p[k], slap_Transpose(work->Qux), work->soln->d[k], -1,
                    1);  // p -= Qux'd
    }
  }
    if (model[0].ltv && !model[0].affine) {
    for (int k = N - 2; k >= 0; --k) {
      // Stage cost expansion
      tiny_ExpandStageCost(work, k);
      // State Gradient: Qx = q + A'p
      slap_MatMulAdd(work->Qx, slap_Transpose(model[0].A[k]), work->soln->p[k + 1], 1, 1);

      // Control Gradient: Qu = r + B'p
      slap_MatMulAdd(work->Qu, slap_Transpose(model[0].B[k]), work->soln->p[k + 1], 1, 1);

      // State Hessian: Qxx = Q + A'P*A
      slap_MatMulAdd(work->soln->P[k], work->soln->P[k + 1], model[0].A[k], 1,
                    0);  // P[k] = P[k+1]*A
      slap_MatMulAdd(work->Qxx, slap_Transpose(model[0].A[k]), work->soln->P[k], 1,
                    1);  // Qxx = Q + A'P*A

      // Control Hessian Quu = R + B'P*B
      slap_MatMulAdd(work->Qxu, work->soln->P[k + 1], model[0].B[k], 1, 0);       // Qxu = P * B
      slap_MatMulAdd(work->Quu, slap_Transpose(model[0].B[k]), work->Qxu, 1, 1);  // Quu = R + B'P*B
      slap_MatMulAdd(work->Quu, slap_Transpose(model[0].B[k]), model[0].B[k], work->reg, 1);
      // Hessian Cross-Term
      slap_MatMulAdd(work->Qux, slap_Transpose(model[0].B[k]), work->soln->P[k], 1,
                    0);  // Qux = B'P*A
      //========= Control constraints ==========
      if (work->stgs->en_cstr_inputs) {
        tiny_EvalInputConstraint(work, k);  // work->cu size = 2*NINPUTS
        tiny_ActiveIneqMask(&(work->cu_mask), work->soln->YU[k], work->cu);
        slap_ScaleByConst(work->cu_mask, work->penalty);  // mask = ρ*mask
        // tiny_EvalInputConstraintJacobian(&cu_jac, *prob);
        // Qu  += G'*(μ[k] + (mask * g)
        // tiny_EvalInputConstraintOffset(&work->cu, *prob);  // g
        slap_Copy(work->cu2, work->soln->YU[k]);
        slap_MatMulAdd(work->cu2, work->cu_mask, work->data->bcu, -1, 1);
        slap_MatMulAdd(work->Qu, slap_Transpose(work->data->Acu), work->cu2, 1, 1);
        // Quu += ∇hu'*mask*∇hu
        slap_MatMulAdd(work->cu_jac2, work->cu_mask, work->data->Acu, 1, 0);
        slap_MatMulAdd(work->Quu, slap_Transpose(work->data->Acu), work->cu_jac2, 1, 1);
      }

      //========= State constraints ==========
      if (work->stgs->en_cstr_states) {
        tiny_EvalStateConstraint(work, k);  // cx size = 2*NINPUTS
        tiny_ActiveIneqMask(&(work->cx_mask), work->soln->YX[k], work->cx);
        slap_ScaleByConst(work->cx_mask, work->penalty);  // mask = ρ*mask
        // tiny_EvalStateConstraintJacobian(&cx_jac, *prob);
        // Qx  += G'*(μx[k] - ρ*mask * g)
        // tiny_EvalStateConstraintOffset(&cx, *prob);  // g
        slap_Copy(work->cx2, work->soln->YX[k]);
        slap_MatMulAdd(work->cx2, work->cx_mask, work->data->bcx, -1, 1);
        slap_MatMulAdd(work->Qx, slap_Transpose(work->data->Acx), work->cx2, 1, 1);
        // Qxx += ρ*∇hx'*mask*∇hx
        slap_MatMulAdd(work->cx_jac2, work->cx_mask, work->data->Acx, 1, 0);
        slap_MatMulAdd(work->Qxx, slap_Transpose(work->data->Acx), work->cx_jac2, 1, 1);
      }
      // Calculate Gains
      slap_Copy(Quu_temp, work->Quu);
      slap_Cholesky(Quu_temp);
      slap_Copy(work->soln->K[k], work->Qux);
      slap_Copy(work->soln->d[k], work->Qu);
      slap_CholeskySolve(Quu_temp, work->soln->d[k]);  // d = Quu\Qu
      slap_CholeskySolve(Quu_temp, work->soln->K[k]);  // K = Quu\Qux

      // Cost-to-Go Hessian: P = Qxx + K'Quu*K - K'Qux - Qux'K
      slap_Copy(work->soln->P[k], work->Qxx);                                  // P = Qxx
      slap_MatMulAdd(work->Qxu, slap_Transpose(work->soln->K[k]), work->Quu, 1, 0);  // Qxu = K'Quu
      slap_MatMulAdd(work->soln->P[k], work->Qxu, work->soln->K[k], 1, 1);           // P += K'Quu*K
      slap_MatMulAdd(work->soln->P[k], slap_Transpose(work->soln->K[k]), work->Qux, -2,
                    1);  // P -= K'Qux
      // slap_MatMulAdd(work->soln->P[k], slap_Transpose(Qux), work->soln->K[k], -1,
      //                1);  // P -= Qux'K

      // Cost-to-Go Gradient: p = Qx + K'Quu*d - K'Qu - Qux'd
      slap_Copy(work->soln->p[k], work->Qx);                          // p = Qx
      slap_MatMulAdd(work->soln->p[k], work->Qxu, work->soln->d[k], 1, 1);  // p += K'Quu*d
      slap_MatMulAdd(work->soln->p[k], slap_Transpose(work->soln->K[k]), work->Qu, -1,
                    1);  // p -= K'Qu
      slap_MatMulAdd(work->soln->p[k], slap_Transpose(work->Qux), work->soln->d[k], -1,
                    1);  // p -= Qux'd
    }
  }
  if (!model[0].ltv && model[0].affine) {
    // for (int k = N - 2; k >= N-2; --k) {
    for (int k = N - 2; k >= 0; --k) {
      // Stage cost expansion
      tiny_ExpandStageCost(work, k);
      // State Gradient: Qx = q + A'(P*f + p)
      slap_Copy(work->soln->p[k], work->soln->p[k + 1]);
      // PrintMatrix(work->soln->p[k]);
      slap_MatMulAdd(work->soln->p[k], work->soln->P[k + 1], model[0].f[0], 1,
                    1);  // p[k] = P[k+1]*f + p[k+1]
      slap_MatMulAdd(work->Qx, slap_Transpose(model[0].A[0]), work->soln->p[k], 1, 1);
      // PrintMatrix(work->Qx);
      // Control Gradient: Qu = r + B'(P*f + p)
      slap_MatMulAdd(work->Qu, slap_Transpose(model[0].B[0]), work->soln->p[k], 1, 1);

      // State Hessian: Qxx = Q + A'P*A
      slap_MatMulAdd(work->soln->P[k], work->soln->P[k + 1], model[0].A[0], 1,
                    0);  // P[k] = P[k+1]*A
      slap_MatMulAdd(work->Qxx, slap_Transpose(model[0].A[0]), work->soln->P[k], 1,
                    1);  // Qxx = Q + A'P*A

      // Control Hessian Quu = R + B'P*B
      slap_MatMulAdd(work->Qxu, work->soln->P[k + 1], model[0].B[0], 1, 0);       // Qxu = P * B
      slap_MatMulAdd(work->Quu, slap_Transpose(model[0].B[0]), work->Qxu, 1, 1);  // Quu = R + B'P*B
      slap_MatMulAdd(work->Quu, slap_Transpose(model[0].B[0]), model[0].B[0], work->reg, 1);
      // Hessian Cross-Term
      slap_MatMulAdd(work->Qux, slap_Transpose(model[0].B[0]), work->soln->P[k], 1,
                    0);  // Qux = B'P*A
      //========= Control constraints ==========
      if (work->stgs->en_cstr_inputs) {
        tiny_EvalInputConstraint(work, k);  // work->cu size = 2*NINPUTS
        tiny_ActiveIneqMask(&(work->cu_mask), work->soln->YU[k], work->cu);
        slap_ScaleByConst(work->cu_mask, work->penalty);  // mask = ρ*mask
        // tiny_EvalInputConstraintJacobian(&cu_jac, *prob);
        // Qu  += G'*(μ[k] + (mask * g)
        // tiny_EvalInputConstraintOffset(&work->cu, *prob);  // g
        slap_Copy(work->cu2, work->soln->YU[k]);
        slap_MatMulAdd(work->cu2, work->cu_mask, work->data->bcu, -1, 1);
        slap_MatMulAdd(work->Qu, slap_Transpose(work->data->Acu), work->cu2, 1, 1);
        // Quu += ∇hu'*mask*∇hu
        slap_MatMulAdd(work->cu_jac2, work->cu_mask, work->data->Acu, 1, 0);
        slap_MatMulAdd(work->Quu, slap_Transpose(work->data->Acu), work->cu_jac2, 1, 1);
      }

      //========= State constraints ==========
      if (work->stgs->en_cstr_states) {
        tiny_EvalStateConstraint(work, k);  // cx size = 2*NINPUTS
        tiny_ActiveIneqMask(&(work->cx_mask), work->soln->YX[k], work->cx);
        slap_ScaleByConst(work->cx_mask, work->penalty);  // mask = ρ*mask
        // tiny_EvalStateConstraintJacobian(&cx_jac, *prob);
        // Qx  += G'*(μx[k] - ρ*mask * g)
        // tiny_EvalStateConstraintOffset(&cx, *prob);  // g
        slap_Copy(work->cx2, work->soln->YX[k]);
        slap_MatMulAdd(work->cx2, work->cx_mask, work->data->bcx, -1, 1);
        slap_MatMulAdd(work->Qx, slap_Transpose(work->data->Acx), work->cx2, 1, 1);
        // Qxx += ρ*∇hx'*mask*∇hx
        slap_MatMulAdd(work->cx_jac2, work->cx_mask, work->data->Acx, 1, 0);
        slap_MatMulAdd(work->Qxx, slap_Transpose(work->data->Acx), work->cx_jac2, 1, 1);
      }
      // Calculate Gains
      slap_Copy(Quu_temp, work->Quu);
      slap_Cholesky(Quu_temp);
      slap_Copy(work->soln->K[k], work->Qux);
      slap_Copy(work->soln->d[k], work->Qu);
      slap_CholeskySolve(Quu_temp, work->soln->d[k]);  // d = Quu\Qu
      slap_CholeskySolve(Quu_temp, work->soln->K[k]);  // K = Quu\Qux

      // Cost-to-Go Hessian: P = Qxx + K'Quu*K - K'Qux - Qux'K
      slap_Copy(work->soln->P[k], work->Qxx);                                  // P = Qxx
      slap_MatMulAdd(work->Qxu, slap_Transpose(work->soln->K[k]), work->Quu, 1, 0);  // Qxu = K'Quu
      slap_MatMulAdd(work->soln->P[k], work->Qxu, work->soln->K[k], 1, 1);           // P += K'Quu*K
      slap_MatMulAdd(work->soln->P[k], slap_Transpose(work->soln->K[k]), work->Qux, -2,
                    1);  // P -= K'Qux
      // slap_MatMulAdd(work->soln->P[k], slap_Transpose(Qux), work->soln->K[k], -1,
      //                1);  // P -= Qux'K

      // Cost-to-Go Gradient: p = Qx + K'Quu*d - K'Qu - Qux'd
      slap_Copy(work->soln->p[k], work->Qx);                          // p = Qx
      slap_MatMulAdd(work->soln->p[k], work->Qxu, work->soln->d[k], 1, 1);  // p += K'Quu*d
      slap_MatMulAdd(work->soln->p[k], slap_Transpose(work->soln->K[k]), work->Qu, -1,
                    1);  // p -= K'Qu
      slap_MatMulAdd(work->soln->p[k], slap_Transpose(work->Qux), work->soln->d[k], -1,
                    1);  // p -= Qux'd
    }
  }
  if (!model[0].ltv && !model[0].affine) {
    for (int k = N - 2; k >= 0; --k) {
      // Stage cost expansion
      tiny_ExpandStageCost(work, k);
      // State Gradient: Qx = q + A'p
      slap_MatMulAdd(work->Qx, slap_Transpose(model[0].A[0]), work->soln->p[k + 1], 1, 1);
      // Control Gradient: Qu = r + B'p
      slap_MatMulAdd(work->Qu, slap_Transpose(model[0].B[0]), work->soln->p[k + 1], 1, 1);

      // State Hessian: Qxx = Q + A'P*A
      slap_MatMulAdd(work->soln->P[k], work->soln->P[k + 1], model[0].A[0], 1,
                    0);  // P[k] = P[k+1]*A
      slap_MatMulAdd(work->Qxx, slap_Transpose(model[0].A[0]), work->soln->P[k], 1,
                    1);  // Qxx = Q + A'P*A

      // Control Hessian Quu = R + B'P*B
      slap_MatMulAdd(work->Qxu, work->soln->P[k + 1], model[0].B[0], 1, 0);       // Qxu = P * B
      slap_MatMulAdd(work->Quu, slap_Transpose(model[0].B[0]), work->Qxu, 1, 1);  // Quu = R + B'P*B
      slap_MatMulAdd(work->Quu, slap_Transpose(model[0].B[0]), model[0].B[0], work->reg, 1);
      // Hessian Cross-Term
      slap_MatMulAdd(work->Qux, slap_Transpose(model[0].B[0]), work->soln->P[k], 1,
                    0);  // Qux = B'P*A
      //========= Control constraints ==========
      if (work->stgs->en_cstr_inputs) {
        tiny_EvalInputConstraint(work, k);  // work->cu size = 2*NINPUTS
        tiny_ActiveIneqMask(&(work->cu_mask), work->soln->YU[k], work->cu);
        slap_ScaleByConst(work->cu_mask, work->penalty);  // mask = ρ*mask
        // tiny_EvalInputConstraintJacobian(&cu_jac, *prob);
        // Qu  += G'*(μ[k] + (mask * g)
        // tiny_EvalInputConstraintOffset(&work->cu, *prob);  // g
        slap_Copy(work->cu2, work->soln->YU[k]);
        slap_MatMulAdd(work->cu2, work->cu_mask, work->data->bcu, -1, 1);
        slap_MatMulAdd(work->Qu, slap_Transpose(work->data->Acu), work->cu2, 1, 1);
        // Quu += ∇hu'*mask*∇hu
        slap_MatMulAdd(work->cu_jac2, work->cu_mask, work->data->Acu, 1, 0);
        slap_MatMulAdd(work->Quu, slap_Transpose(work->data->Acu), work->cu_jac2, 1, 1);
      }

      //========= State constraints ==========
      if (work->stgs->en_cstr_states) {
        tiny_EvalStateConstraint(work, k);  // cx size = 2*NINPUTS
        tiny_ActiveIneqMask(&(work->cx_mask), work->soln->YX[k], work->cx);
        slap_ScaleByConst(work->cx_mask, work->penalty);  // mask = ρ*mask
        // tiny_EvalStateConstraintJacobian(&cx_jac, *prob);
        // Qx  += G'*(μx[k] - ρ*mask * g)
        // tiny_EvalStateConstraintOffset(&cx, *prob);  // g
        slap_Copy(work->cx2, work->soln->YX[k]);
        slap_MatMulAdd(work->cx2, work->cx_mask, work->data->bcx, -1, 1);
        slap_MatMulAdd(work->Qx, slap_Transpose(work->data->Acx), work->cx2, 1, 1);
        // Qxx += ρ*∇hx'*mask*∇hx
        slap_MatMulAdd(work->cx_jac2, work->cx_mask, work->data->Acx, 1, 0);
        slap_MatMulAdd(work->Qxx, slap_Transpose(work->data->Acx), work->cx_jac2, 1, 1);
      }
      // Calculate Gains
      slap_Copy(Quu_temp, work->Quu);
      slap_Cholesky(Quu_temp);
      slap_Copy(work->soln->K[k], work->Qux);
      slap_Copy(work->soln->d[k], work->Qu);
      slap_CholeskySolve(Quu_temp, work->soln->d[k]);  // d = Quu\Qu
      slap_CholeskySolve(Quu_temp, work->soln->K[k]);  // K = Quu\Qux

      // Cost-to-Go Hessian: P = Qxx + K'Quu*K - K'Qux - Qux'K
      slap_Copy(work->soln->P[k], work->Qxx);                                  // P = Qxx
      slap_MatMulAdd(work->Qxu, slap_Transpose(work->soln->K[k]), work->Quu, 1, 0);  // Qxu = K'Quu
      slap_MatMulAdd(work->soln->P[k], work->Qxu, work->soln->K[k], 1, 1);           // P += K'Quu*K
      slap_MatMulAdd(work->soln->P[k], slap_Transpose(work->soln->K[k]), work->Qux, -2,
                    1);  // P -= K'Qux
      // slap_MatMulAdd(work->soln->P[k], slap_Transpose(Qux), work->soln->K[k], -1,
      //                1);  // P -= Qux'K

      // Cost-to-Go Gradient: p = Qx + K'Quu*d - K'Qu - Qux'd
      slap_Copy(work->soln->p[k], work->Qx);                          // p = Qx
      slap_MatMulAdd(work->soln->p[k], work->Qxu, work->soln->d[k], 1, 1);  // p += K'Quu*d
      slap_MatMulAdd(work->soln->p[k], slap_Transpose(work->soln->K[k]), work->Qu, -1,
                    1);  // p -= K'Qu
      slap_MatMulAdd(work->soln->p[k], slap_Transpose(work->Qux), work->soln->d[k], -1,
                    1);  // p -= Qux'd
    }
  }  
  // tiny_ExpandTerminalCost(work);
  slap_Copy(work->soln->P[N-1], work->data->Qf);
  // Replace P[N] since we used it for Quu_temp (need improving later)
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SolveAlLqr(tiny_Workspace* work) {
  slap_Copy(work->soln->X[0], work->data->x0);
  // Shortcut unconstrained problem
  if (!IsConstrained(work)) {
    printf("UNCONSTRAINED!\n");
    tiny_ForwardPass(work);
    tiny_BackwardPass(work);
    return TINY_NO_ERROR;
  }

  // int exitflag;
  int verbose = work->stgs->verbose;         // Boolean whether you can print
  // int compute_cost_function; // Boolean: compute the cost function in the loop or not

  tiny_Model* model = work->data->model;
  int N = model[0].nhorizon;
  // int n = model[0].nstates;
  // int m = model[0].ninputs;
 
  tiny_ResetInfo(work);

  // Roll-out initial guess
  tiny_RollOutOpenLoop(work);

  for (work->info->iter_al = 1; work->info->iter_al < work->stgs->max_iter_al;
       ++work->info->iter_al) {
    for (int iter = 1; iter < work->stgs->max_iter_riccati; ++iter) {
      work->info->iter_riccati += 1;
      if (verbose > 1) printf("backward pass\n");
      tiny_ConstrainedBackwardPass(work);
      if (verbose > 1) printf("forward pass\n");
      tiny_ConstrainedForwardPass(work);
    }

    sfloat norm_inf = 0.0;   // temporary var to save inf norm of vector

    if (work->stgs->en_cstr_inputs) {
      for (int k = 0; k < N - 1; ++k) {
        //========= Control constraints ==========
        tiny_EvalInputConstraint(work, k);  // work->cu size = 2*NINPUTS
        tiny_ActiveIneqMask(&(work->cu_mask), work->soln->YU[k], work->cu);
        slap_ScaleByConst(work->cu_mask, work->penalty);  // mask = ρ*mask
        // Constraint violation
        slap_ArgMax(work->cu, &norm_inf);
        norm_inf = norm_inf > 0.0 ? norm_inf : 0.0;
        norm_inf = norm_inf * 2;
        // convio = max(convio,norm(hxv + abs.(hxv),Inf))
        work->info->dua_res = work->info->dua_res < norm_inf ? norm_inf : work->info->dua_res;
        // Update duals
        slap_Copy(work->YU_hat, work->soln->YU[k]);
        slap_MatMulAdd(work->YU_hat, work->cu_mask, work->data->bcu, -1, 1);  //μ[k] - ρ*mask * g
        tiny_ProjectOrthantDuals(&(work->soln->YU[k]), work->YU_hat);
      }
    }

    if (work->stgs->en_cstr_states) {
      for (int k = 0; k < N; ++k) {
        //========= State constraints ==========
        tiny_EvalStateConstraint(work, k);  // work->cu size = 2*NINPUTS
        tiny_ActiveIneqMask(&(work->cx_mask), work->soln->YX[k], work->cx);
        slap_ScaleByConst(work->cx_mask, work->penalty);  // mask = ρ*mask
        // Constraint violation
        slap_ArgMax(work->cx, &norm_inf);
        norm_inf = norm_inf > 0.0 ? norm_inf : 0.0;
        norm_inf = norm_inf * 2;
        // convio = max(convio,norm(hxv + abs.(hxv),Inf))
        work->info->dua_res = work->info->dua_res < norm_inf ? norm_inf : work->info->dua_res;
        // Update duals
        // tiny_IneqStatesOffset(&ineq_state, *prob);  // g
        slap_Copy(work->YX_hat, work->soln->YX[k]);
        slap_MatMulAdd(work->YX_hat, work->cx_mask, work->data->bcx, -1, 1);  //μ[k] - ρ*mask * g
        tiny_ProjectOrthantDuals(&(work->soln->YX[k]), work->YX_hat);
      }
    }

    if (work->stgs->en_cstr_goal) {
      //========= Goal constraints ==========
      slap_MatrixAddition(work->cg, work->soln->X[N - 1], work->data->X_ref[N - 1], -1);
      norm_inf = slap_NormInf(work->cg);
      work->info->dua_res = work->info->dua_res < norm_inf ? norm_inf : work->info->dua_res;
      // λ -= ρ*h
      slap_Copy(work->cg, work->data->X_ref[N - 1]);  // h = xg
      slap_MatrixAddition(work->soln->YG, work->soln->YG, work->cg, -work->penalty);
    }    

    // printing
    if (verbose > 0) {
      if (work->info->iter_al == 1) {
        printf("iter           J           d      convio       reg         ρ\n");
        printf("------------------------------------------------------------\n");
      }
      printf("%4d%12.3e%12.3e%12.3e%10.1e%10.1e\n", 
             work->info->iter_al, work->info->obj_al, work->info->pri_res,
             work->info->dua_res, work->reg, work->penalty);
    }

    // check AL termination
    if (work->info->dua_res < work->stgs->tol_abs_cstr) {
      if (verbose > 0) printf("SUCCESS!\n");
      work->penalty = work->stgs->penalty_init;  // reset penalty for next MPC
      return TINY_NO_ERROR;
    }

    // update penalty
    work->penalty = work->penalty * work->stgs->penalty_mul;    
  }
  work->info->status_val = MAX_ITER_AL;
  return TINY_NO_ERROR;
}