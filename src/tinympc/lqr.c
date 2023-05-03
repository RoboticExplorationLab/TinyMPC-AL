#include "lqr.h"

enum tiny_ErrorCode tiny_RollOutClosedLoopCost(tiny_Workspace* work) {
  tiny_Model* model = work->data->model;
  int N = model[0].nhorizon;
  int adaptive_horizon = work->stgs->adaptive_horizon;
  work->info->obj_pri = 0.0;
  
  if (model[0].ltv) {
    for (int k = 0; k < N - 1; ++k) {
      // Control input: u = - d - K*x
      slap_Copy(work->soln->U[k], work->soln->d[k]); // u[k] = -d[k]
      slap_MatMulAdd(work->soln->U[k], work->soln->K[k], work->soln->X[k], -1, -1);  // u[k] -= K[k] * x[k]
      // Next state: x = A*x + B*u + f
      if (adaptive_horizon && k > adaptive_horizon - 1) {
        tiny_EvalModel(&(work->soln->X[k + 1]), work->soln->X[k], work->soln->U[k], &model[1], k);
      }
      else {
        tiny_EvalModel(&(work->soln->X[k + 1]), work->soln->X[k], work->soln->U[k], &model[0], k);
      }
      tiny_AddStageCost(work, k);
    }    
    tiny_AddTerminalCost(work);
  }
  else {
    for (int k = 0; k < N - 1; ++k) {
      // Control input: u = - d - K*x
      slap_Copy(work->soln->U[k], work->soln->d[k]); // u[k] = -d[k]
      slap_MatMulAdd(work->soln->U[k], work->soln->K[k], work->soln->X[k], -1, -1);  // u[k] -= K[k] * x[k]
      // Next state: x = A*x + B*u + f
      if (adaptive_horizon && k > adaptive_horizon - 1) {        
        tiny_EvalModel(&(work->soln->X[k + 1]), work->soln->X[k], work->soln->U[k], &model[1], 0);
      }
      else {
        tiny_EvalModel(&(work->soln->X[k + 1]), work->soln->X[k], work->soln->U[k], &model[0], 0);
      }
      tiny_AddStageCost(work, k);
      // printf("%f\n", work->info->obj_pri);
    } 
    // printf("%f\n", work->info->obj_pri);
    tiny_AddTerminalCost(work);       
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_ForwardPass(tiny_Workspace* work) {
  // tiny_RollOutClosedLoop(work);
  tiny_RollOutClosedLoopCost(work);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_BackwardPass(tiny_Workspace* work) {
  tiny_Model* model = work->data->model;
  int N = model[0].nhorizon;
  // int n = model[0].nstates;
  int m = model[0].ninputs;

  tiny_ExpandTerminalCost(work);
  // Hijack the first part of P[N] to use for Cholesky decomposition
  // NOTE: Assumes m <= n
  Matrix Quu_temp = slap_Reshape(work->soln->P[N - 1], m, m);

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
    // printf("LTI AND AFFINE\n");
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


enum tiny_ErrorCode tiny_SolveLqr(tiny_Workspace* work) {

  slap_Copy(work->soln->X[0], work->data->x0);

  tiny_ResetInfo(work);

  tiny_BackwardPass(work);
  tiny_ForwardPass(work);

  work->info->status_val = TINY_SOLVED;
  
  return TINY_NO_ERROR;
}