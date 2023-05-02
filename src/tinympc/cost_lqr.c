#include "cost_lqr.h"

enum tiny_ErrorCode tiny_AddStageCost(tiny_Workspace* work, const int k) {
  int n = work->data->model[0].nstates;
  int m = work->data->model[0].ninputs;
  sfloat dx_data[n];
  Matrix dx = slap_MatrixFromArray(n, 1, dx_data);
  slap_MatrixAddition(dx, work->soln->X[k], work->data->X_ref[k], -1);
  work->info->obj_pri += 0.5 * slap_QuadraticForm(dx, work->data->Q, dx);
  // PrintMatrixT(dx);
  // printf("%f\n", work->info->obj_pri);
  Matrix du = slap_MatrixFromArray(m, 1, dx_data);
  slap_MatrixAddition(du, work->soln->U[k], work->data->U_ref[k], -1);
  work->info->obj_pri += 0.5 * slap_QuadraticForm(du, work->data->R, du);

  // PrintMatrixT(du);
  // PrintMatrix(work->data->R);
  // printf("%f\n", work->info->obj_pri);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_AddTerminalCost(tiny_Workspace* work) {
  int n = work->data->model[0].nstates;
  int N = work->data->model[0].nhorizon;
  sfloat dx_data[n];
  Matrix dx = slap_MatrixFromArray(n, 1, dx_data);
  slap_MatrixAddition(dx, work->soln->X[N - 1], work->data->X_ref[N - 1], -1);
  work->info->obj_pri += 0.5 * slap_QuadraticForm(dx, work->data->Qf, dx);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_ExpandStageCost(tiny_Workspace* work, const int k) {
  slap_Copy(work->Qxx, work->data->Q);
  slap_Copy(work->Quu, work->data->R);

  slap_Copy(work->Qx, work->data->q[k]);
  slap_Copy(work->Qu, work->data->r[k]); 
  return TINY_NO_ERROR; 
}

enum tiny_ErrorCode tiny_ExpandTerminalCost(tiny_Workspace* work) {
  int N = work->data->model[0].nhorizon;
  slap_Copy(work->soln->P[N-1], work->data->Qf);
  slap_Copy(work->soln->p[N-1], work->data->qf);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_UpdateLinearCost(tiny_Workspace* work) {
  int N = work->data->model[0].nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    slap_MatMulAdd(work->data->q[k], work->data->Q, work->data->X_ref[k], -1, 0);
    slap_MatMulAdd(work->data->r[k], work->data->R, work->data->U_ref[k], -1, 0);
  }
  slap_MatMulAdd(work->data->qf, work->data->Qf, work->data->X_ref[N-1], -1, 0);
  return TINY_NO_ERROR;
}