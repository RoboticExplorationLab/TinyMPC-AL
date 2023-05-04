#include "auxil.h"

enum tiny_ErrorCode tiny_InitSettings(tiny_Settings* stgs) {
  SLAP_ASSERT(stgs != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_InitSettings: settings must not be TINY_NULL");
  stgs->reg_min       = (sfloat)REG_MIN;
  stgs->reg_max       = (sfloat)REG_MAX;
  stgs->reg_mul       = (sfloat)REG_MUL;
  stgs->en_reg_update = EN_REG_UPDATE;

  stgs->penalty_init  = (sfloat)PENALTY_INIT;
  stgs->penalty_max   = (sfloat)PENALTY_MAX;
  stgs->penalty_mul   = (sfloat)PENALTY_MUL;

  stgs->alpha_mul     = (sfloat)ALPHA_MUL;

  stgs->max_iter_al   = MAX_ITER_AL;
  stgs->max_iter_riccati = MAX_ITER_RICCATI;
  stgs->max_iter_ls   = MAX_ITER_LS;

  stgs->tol_abs_riccati = (sfloat)TOL_ABS_RICCATI;
  stgs->tol_abs_cstr    = (sfloat)TOL_ABS_CSTR;

  stgs->en_cstr_states  = EN_CSTR_STATES;
  stgs->en_cstr_inputs  = EN_CSTR_INPUTS;
  stgs->en_cstr_goal    = EN_CSTR_GOAL;

  stgs->verbose          = VERBOSE;
  stgs->adaptive_horizon = ADAPTIVE_HORIZON;
  stgs->check_riccati    = CHECK_RICCATI;
  stgs->check_al         = CHECK_AL;
  stgs->warm_start       = WARM_START;
  stgs->time_limit       = TIME_LIMIT;

  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetUnconstrained(tiny_Settings* stgs) {
  stgs->en_cstr_states  = 0;
  stgs->en_cstr_inputs  = 0;
  stgs->en_cstr_goal    = 0;
  stgs->check_al        = 0;

  return TINY_NO_ERROR;
}

// enum tiny_ErrorCode tiny_InitSolution(tiny_Workspace* work) {
//   tiny_Solution* soln = work->soln;
//   tiny_Model* model   = &(work->data->model[0]);
//   int n = model->nstates;
//   int m = model->ninputs;
//   int N = model->nhorizon;

//   soln->X = TINY_NULL;
//   soln->U = TINY_NULL;

//   soln->YX = TINY_NULL;
//   soln->YU = TINY_NULL;
//   soln->YG = TINY_NULL_MAT;

//   soln->K = TINY_NULL;
//   soln->d = TINY_NULL;
//   soln->P = TINY_NULL;
//   soln->p = TINY_NULL;

//   soln->data_size = n*N + m*(N-1) + (m+1)*n*(N-1) + (n+1)*n*N;

//   if (work->stgs->en_cstr_inputs) {
//     soln->data_size += m*(N-1)*2;
//   }

//   if (work->stgs->en_cstr_states) {
//     soln->data_size += n*N*2;
//   }  

//   if (work->stgs->en_cstr_goal) {
//     soln->data_size += n;
//   }  
//   return TINY_NO_ERROR;
// }

enum tiny_ErrorCode tiny_InitSolutionFromMatrix(tiny_Workspace* work, 
                                               Matrix* X, Matrix* U,
                                               Matrix* P, Matrix* p,
                                               Matrix* K, Matrix* d,
                                               Matrix* YX, Matrix* YU,
                                               Matrix YG) {
  work->soln->X = X;
  work->soln->U = U;

  work->soln->P = P;
  work->soln->p = p;
  work->soln->K = K;
  work->soln->d = d;

  // if (work->stgs->en_cstr_states) {
  //   SLAP_ASSERT(YX != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  //   "tiny_InitSolitionDataArray: YX must not be TINY_NULL");
  // }
  work->soln->YX = YX;

  // if (work->stgs->en_cstr_inputs) {
  //   SLAP_ASSERT(YU != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  //   "tiny_InitSolitionDataArray: YU must not be TINY_NULL");
  // } 
  work->soln->YU = YU;

  // if (work->stgs->en_cstr_goal) {
  //   SLAP_ASSERT(slap_IsNull(YG), SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  //   "tiny_InitSolitionDataArray: YG must not be TINY_NULL");
  // }
  work->soln->YG = YG;

  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitSolnTrajFromArray(tiny_Workspace* work,
Matrix* X, Matrix* U,
sfloat* X_data, sfloat* U_data) {

  int N = work->data->model->nhorizon;
  int n = work->data->model->nstates;
  int m = work->data->model->ninputs;

  work->soln->X = X;  
  work->soln->U = U;

  for (int i = 0; i < N - 1; ++i) {
    U[i] = slap_MatrixFromArray(m, 1, &U_data[i * m]);
  }

  for (int i = 0; i < N; ++i) {
    X[i] = slap_MatrixFromArray(n, 1, &X_data[i * n]);
  }

  return TINY_NO_ERROR;  
}

enum tiny_ErrorCode tiny_InitSolnDualsFromArray(tiny_Workspace* work,
Matrix* YX, Matrix* YU,
sfloat* YX_data, sfloat* YU_data, sfloat* YG_data) {

  int N = work->data->model->nhorizon;
  int n = work->data->model->nstates;
  int m = work->data->model->ninputs;

  work->soln->YX = YX;  
  work->soln->YU = YU;
  if (YG_data) {
    work->soln->YG = slap_MatrixFromArray(n, 1, YG_data);
  }

  if (YU_data) {
    for (int i = 0; i < N - 1; ++i) {
      YU[i] = slap_MatrixFromArray(2 * m, 1, &YU_data[i * 2 * m]);
    }
  }

  if (YX_data) {
    for (int i = 0; i < N; ++i) {
      YX[i] = slap_MatrixFromArray(2 * n, 1, &YX_data[2 * n]);
    }
  }

  return TINY_NO_ERROR;  
}

enum tiny_ErrorCode tiny_InitSolnGainsFromArray(tiny_Workspace* work, 
Matrix* K, Matrix* d, Matrix* P, Matrix* p, 
sfloat* K_data, sfloat* d_data, sfloat* P_data, sfloat* p_data) {

  int N = work->data->model->nhorizon;
  int n = work->data->model->nstates;
  int m = work->data->model->ninputs;

  work->soln->K = K;  
  work->soln->d = d;
  work->soln->P = P;  
  work->soln->p = p;

  for (int i = 0; i < N - 1; ++i) {
    K[i] = slap_MatrixFromArray(m, n, &K_data[i * m * n]);
    d[i] = slap_MatrixFromArray(m, 1, &d_data[i * m]);
  }
  for (int i = 0; i < N; ++i) {
    P[i] = slap_MatrixFromArray(n, n, &P_data[i * n * n]);
    p[i] = slap_MatrixFromArray(n, 1, &p_data[i * n]);
  }

  return TINY_NO_ERROR;  
}

// enum tiny_ErrorCode tiny_InitData(tiny_Workspace* work) {
//   tiny_Data* data   = work->data;
//   tiny_Model* model = &(work->data->model[0]);
//   int n = model->nstates;
//   int m = model->ninputs;
//   int N = model->nhorizon;

//   data->x0 = TINY_NULL_MAT;
//   data->Q = TINY_NULL_MAT;
//   data->R = TINY_NULL_MAT;
//   data->Qf = TINY_NULL_MAT;
//   data->q = TINY_NULL;
//   data->r = TINY_NULL;
//   data->qf = TINY_NULL_MAT;
//   data->X_ref = TINY_NULL;
//   data->U_ref = TINY_NULL;
//   data->Acx = TINY_NULL_MAT;
//   data->bcx = TINY_NULL_MAT;
//   data->Acu = TINY_NULL_MAT;
//   data->bcu = TINY_NULL_MAT;

//   data->data_size = n + n*n*2 + m*m + (N-1)*(n + m) + m + N*n + (N-1)*m;  // with ref

//   if (work->stgs->en_cstr_inputs) {
//     data->data_size += 2*m*m + 2*m; 
//   }

//   if (work->stgs->en_cstr_states) {
//     data->data_size += 2*n*n + 2*n;
//   }  

//   return TINY_NO_ERROR;
// }

enum tiny_ErrorCode tiny_InitDataFromMatrix(tiny_Workspace* work, Matrix x0,
                                            Matrix Q, Matrix R, Matrix Qf,
                                            Matrix* q, Matrix* r, Matrix qf,
                                            Matrix* X_ref, Matrix* U_ref,
                                            Matrix Acx, Matrix bcx,
                                            Matrix Acu, Matrix bcu) {
  work->data->x0 = x0;

  work->data->Q = Q;
  work->data->R = R;
  work->data->Qf = Qf;
  work->data->q = q;
  work->data->r = r;
  work->data->qf = qf;

  work->data->X_ref = X_ref;
  work->data->U_ref = U_ref;

  // if (work->stgs->en_cstr_states) {
  //   SLAP_ASSERT(slap_IsNull(Acx) || slap_IsNull(bcx), SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  //   "tiny_InitSolitionDataArray: Acx and bcx must not be TINY_NULL_MAT");
  // }
  work->data->Acx = Acx;
  work->data->bcx = bcx;

  // if (work->stgs->en_cstr_inputs) {
  //   SLAP_ASSERT(slap_IsNull(Acu) || slap_IsNull(bcu), SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  //   "tiny_InitSolitionDataArray: Acu and bcu must not be TINY_NULL_MAT");
  // } 
  work->data->Acu = Acu;
  work->data->bcu = bcu;

  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitDataQuadCostFromArray(tiny_Workspace* work, 
sfloat* Q_data, 
sfloat* R_data, 
sfloat* Qf_data) {

  int n = work->data->model->nstates;
  int m = work->data->model->ninputs;
  work->data->Q = slap_MatrixFromArray(n, n, Q_data);
  work->data->R = slap_MatrixFromArray(m, m, R_data);
  work->data->Qf = slap_MatrixFromArray(n, n, Qf_data);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitDataLinearCostFromArray(tiny_Workspace* work, 
Matrix* q, Matrix* r, 
sfloat* q_data, sfloat* r_data, sfloat* qf_data) {

  int N = work->data->model->nhorizon;
  int n = work->data->model->nstates;
  int m = work->data->model->ninputs;
  work->data->qf = slap_MatrixFromArray(n, 1, qf_data);  
  work->data->q = q;
  work->data->r = r;
  for (int i = 0; i < N - 1; ++i) {
    q[i] = slap_MatrixFromArray(n, 1, &q_data[i * n]);
    r[i] = slap_MatrixFromArray(m, 1, &r_data[i * m]);
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitWorkspace(tiny_Workspace* work,
                                       tiny_Info* info,
                                       tiny_Model* model,
                                       tiny_Data* data,
                                       tiny_Solution* soln,
                                       tiny_Settings* stgs) {
  SLAP_ASSERT(work != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_InitWorkspace: work must not be TINY_NULL");
  SLAP_ASSERT(info != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_InitWorkspace: info must not be TINY_NULL");
  SLAP_ASSERT(model != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_InitWorkspace: model must not be TINY_NULL");
  SLAP_ASSERT(data != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_InitWorkspace: data must not be TINY_NULL");    
  SLAP_ASSERT(stgs != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_InitWorkspace: stgs must not be TINY_NULL");
  SLAP_ASSERT(soln != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_InitWorkspace: soln must not be TINY_NULL");  

  work->data = data;
  work->info = info;
  work->soln = soln;
  work->stgs = stgs;
  work->data->model = model;

  // tiny_InitSolution(work);
  // tiny_InitData(work);

  int n = model->nstates;
  int m = model->ninputs;
  // int N = model->nhorizon;

  work->reg = (sfloat)REG_MIN;
  work->alpha = (sfloat)ALPHA;
  work->penalty = work->stgs->penalty_init;
  // work->Q_temp = TINY_NULL_MAT;
  // work->c_temp = TINY_NULL_MAT;

  // work->Qxx = TINY_NULL_MAT;
  // work->Qxu = TINY_NULL_MAT;
  // work->Qux = TINY_NULL_MAT;
  // work->Quu = TINY_NULL_MAT;
  // work->Qx = TINY_NULL_MAT;
  // work->Qu = TINY_NULL_MAT;

  // work->cu = TINY_NULL_MAT;
  // work->cu2 = TINY_NULL_MAT;
  // work->cu_jac = TINY_NULL_MAT;
  // work->cu_jac2 = TINY_NULL_MAT;
  // work->cu_mask = TINY_NULL_MAT;
  // work->YU_hat = TINY_NULL_MAT;

  // work->cx = TINY_NULL_MAT;
  // work->cx2 = TINY_NULL_MAT;
  // work->cx_jac = TINY_NULL_MAT;
  // work->cx_jac2 = TINY_NULL_MAT;
  // work->cx_mask = TINY_NULL_MAT;
  // work->YX_hat = TINY_NULL_MAT;

  // work->cg = TINY_NULL_MAT;


  work->data_size = 2 * n * (2 * n + 2 * n + 2) + (n + m) * (n + m + 1);
  work->first_run = 1;
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitWorkspaceTempData(tiny_Workspace* work, sfloat* temp_data) {
  int n = work->data->model->nstates;
  int m = work->data->model->ninputs;

  work->Q_temp = slap_MatrixFromArray(n + m, n + m + 1, temp_data);
  work->Qxx = slap_CreateSubMatrix(work->Q_temp, 0, 0, n, n);
  work->Qxu = slap_CreateSubMatrix(work->Q_temp, 0, n, n, m);
  work->Qux = slap_CreateSubMatrix(work->Q_temp, n, 0, m, n);
  work->Quu = slap_CreateSubMatrix(work->Q_temp, n, n, m, m);
  work->Qx = slap_CreateSubMatrix(work->Q_temp, 0, n + m, n, 1);
  work->Qu = slap_CreateSubMatrix(work->Q_temp, n, n + m, m, 1);

  work->c_temp = slap_MatrixFromArray(2 * n, 2 * n + 2 * n + 2,
                                          &temp_data[(n + m) * (n + m + 1)]);
  work->cu = slap_CreateSubMatrix(work->c_temp, 0, 0, 2 * m, 1);
  work->YU_hat = slap_CreateSubMatrix(work->c_temp, 0, 1, 2 * m, 1);
  // work->cu2 = slap_CreateSubMatrix(work->c_temp, 0, 1, 2 * m, 1);   // YU_hat ~ cu2
  work->cu2 = work->YU_hat;
  work->cu_mask = slap_CreateSubMatrix(work->c_temp, 0, 2, 2 * m, 2 * m);
  work->cu_jac = slap_CreateSubMatrix(work->c_temp, 0, 2 * m + 2, 2 * m, m);
  work->cu_jac2 = slap_CreateSubMatrix(work->c_temp, 0, 2 * m + 2 + m, 2 * m, m);

  work->cx = slap_CreateSubMatrix(work->c_temp, 0, 0, 2 * n, 1);
  work->YX_hat = slap_CreateSubMatrix(work->c_temp, 0, 1, 2 * n, 1);  // YX_hat ~ cx2
  // work->cx2 = slap_CreateSubMatrix(work->c_temp, 0, 1, 2 * n, 1);
  work->cx2 = work->YX_hat;
  work->cx_mask = slap_CreateSubMatrix(work->c_temp, 0, 2, 2 * n, 2 * n);
  work->cx_jac = slap_CreateSubMatrix(work->c_temp, 0, 2 * n + 2, 2 * n, n);
  work->cx_jac2 = slap_CreateSubMatrix(work->c_temp, 0, 2 * n + 2 + n, 2 * n, n);

  work->cg = slap_MatrixFromArray(n, 1, &temp_data[2 * n * (2 * n + 2 * n + 2)]);   

  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_ResetInfo(tiny_Workspace* work) {
  work->info->iter_al = 0;
  work->info->iter_riccati = 0;
  work->info->status_val = TINY_UNSOLVED;

  work->info->obj_pri = 0.0;
  work->info->obj_al = 0.0;
  work->info->pri_res = 0.0;
  work->info->dua_res = 0.0;
  return TINY_NO_ERROR;
}