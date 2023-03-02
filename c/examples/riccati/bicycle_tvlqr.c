// Check README.md
// Second-order bicycle model tracking TVLQR in MPC style
// New: Code generation for Jacobians and dynamics
// Task: LTV system to track a reference trajectory
// TODO: convert this to test Riccati

#include <stdio.h>

#include "bicycle.h"
#include "riccati.h"
#include "slap/slap.h"
#include "util.h"

#define H (0.1)
#define NSTATES (5)
#define NINPUTS (2)
#define NHORIZON (10)
#define NSIM (101 + 1 * NHORIZON)  // extend the end to complete MPC

int main(void) {
  printf("\n*** PROBLEM DEFINITION ***\n");
  struct tiny_Model_Bicycle model = tiny_DefaultModel_Bicycle;

  // array data to construct matrix, each column
  double Q_data[NSTATES * NSTATES] = {0};
  double R_data[NINPUTS * NINPUTS] = {0};
  double q_data[NSTATES] = {0};
  double r_data[NINPUTS] = {0};

  // A and B must be lists of matrices
  double A_data[NSTATES * NSTATES] = {0};
  double B_data[NSTATES * NINPUTS] = {0};

  double x0_data[NSTATES] = {1, -1, 0, 1, 1};
  double xf_data[NSTATES] = {20, 20, M_PI / 4, 0, 0};
  double uf_data[NINPUTS] = {0};

  double Pp_data[NSTATES * (NSTATES + 1) * NHORIZON] = {0};  // stores P and p
  double Kd_data[NINPUTS * (NSTATES + 1) * (NHORIZON - 1)] = {
      0};  // stores K and d

  double x_data[NSTATES * NSIM] = {0};
  double u_data[NINPUTS * (NSIM - 1)] = {0};
  double xref_data[NSTATES * NSIM] = {0};
  double uref_data[NINPUTS * (NSIM - 1)] = {0};

  // Read reference trajectory from files
  const char *file_x_ref = "../examples/riccati/data/xref_data.txt";
  const char *file_u_ref = "../examples/riccati/data/uref_data.txt";
  // tiny_ReadData(file_xref, xref_data, NSTATES * NSIM, false);
  // tiny_ReadData(file_uref, uref_data, NINPUTS * (NSIM - 1), false);
  tiny_ReadData_ExtendGoal(file_x_ref, xref_data, xf_data, NSTATES,
                           NSTATES * NSIM, false);
  tiny_ReadData_ExtendGoal(file_u_ref, uref_data, uf_data, NINPUTS,
                           NINPUTS * (NSIM - 1), false);
  // tiny_ReadData_Extend(file_xref, xref_data, NSTATES, NSTATES * NSIM, false);
  // tiny_ReadData_Extend(file_uref, uref_data, NINPUTS, NINPUTS * (NSIM-1),
  // false); Create matrix from array data
  Matrix Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(Q, 1);
  // printf("\nQ = \n"); slap_PrintMatrix(Q);
  tiny_Print(Q);
  Matrix R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(R, 0.01);
  printf("\nR = \n");
  slap_PrintMatrix(R);
  Matrix q = slap_MatrixFromArray(NSTATES, 1, q_data);
  // printf("\nq = \n"); slap_PrintMatrix(q);
  Matrix r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  // printf("\nr = \n"); slap_PrintMatrix(r);
  Matrix A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  Matrix B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);

  Matrix x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  printf("\nx0 = \n");
  slap_PrintMatrix(x0);
  Matrix xf = slap_MatrixFromArray(NSTATES, 1, xf_data);
  tiny_Print(xf);

  // Matrix of pointers
  Matrix P_hist[NHORIZON];
  Matrix p_hist[NHORIZON];
  Matrix K_hist[NHORIZON - 1];
  Matrix d_hist[NHORIZON - 1];
  Matrix x_hist[NSIM];
  Matrix u_hist[NSIM - 1];
  Matrix u_ref[NSIM - 1];
  Matrix x_ref[NSIM];
  // Pointer to the pre-allocated array
  double *Pp = Pp_data;
  double *Kd = Kd_data;
  double *xp = x_data;
  double *up = u_data;
  double *xrefp = xref_data;
  double *urefp = uref_data;

  for (int k = 0; k < NSIM; ++k) {
    // Pointer to each block, then next
    if (k < NHORIZON) {
      P_hist[k] = slap_MatrixFromArray(NSTATES, NSTATES, Pp);
      Pp += NSTATES * NSTATES;
      p_hist[k] = slap_MatrixFromArray(NSTATES, 1, Pp);
      Pp += NSTATES;
      if (k < NHORIZON - 1) {
        K_hist[k] = slap_MatrixFromArray(NINPUTS, NSTATES, Kd);
        Kd += NINPUTS * NSTATES;
        d_hist[k] = slap_MatrixFromArray(NINPUTS, 1, Kd);
        Kd += NINPUTS;
      }
    }
    x_hist[k] = slap_MatrixFromArray(NSTATES, 1, xp);
    xp += NSTATES;
    x_ref[k] = slap_MatrixFromArray(NSTATES, 1, xrefp);
    xrefp += NSTATES;

    if (k < NSIM - 1) {
      u_hist[k] = slap_MatrixFromArray(NINPUTS, 1, up);
      up += NINPUTS;
      u_ref[k] = slap_MatrixFromArray(NINPUTS, 1, urefp);
      urefp += NINPUTS;
    }
  }
  // End of initializing memory and variables

  // Temporary matrix for underlying calculation
  Matrix S = slap_NewMatrixZeros(NSTATES + NINPUTS, NSTATES + NINPUTS + 1);
  slap_MatrixCopy(x_hist[0], x0);
  printf("\n*** START SOLVING ***\n");
  // MPC loop
  for (int k = 0; k < NSIM - NHORIZON; ++k) {
    tiny_Riccati_LTVf(NHORIZON - 1, A, B, tiny_GetJacobians, Q, R, q, r, K_hist,
                      d_hist, P_hist, p_hist, &x_ref[k], &u_ref[k], S);
    // tiny_Print(K_hist[0]);
    // Control input: u = uf - d - K*(x - xf)
    slap_MatrixAddition(u_hist[k], u_ref[k], d_hist[0],
                        -1);  // u[k] = un[k] - d[k]
    // slap_PrintMatrix(slap_Transpose(u_hist[k]));
    slap_MatMulAdd(u_hist[k], K_hist[0], x_hist[k], -1, 1);  // u[k] -= K[k] * x[k]
    // slap_PrintMatrix(slap_Transpose(u_hist[k]));
    slap_MatMulAdd(u_hist[k], K_hist[0], x_ref[k], 1, 1);  // u[k] += K[k] * xn[k]
    // printf("u[%d] = ", k);
    // slap_PrintMatrix(slap_Transpose(u_hist[k]));

    // Next state: x = f(x, u)
    // tiny_Clamps(u_hist[k].data, model.u_min, model.u_max, NINPUTS);
    tiny_Dynamics_RK4_Raw(x_hist[k + 1].data, x_hist[k].data, u_hist[k].data);

    // Tracking errors (show different metrics)
    printf("ex[%d] = %.4f\n", k, slap_MatrixNormedDifference(x_ref[k], x_hist[k]));
    // printf("x[%d] = ", k);
    // slap_PrintMatrix(slap_Transpose(u_hist[k]));
  }
  printf("\n*** END OF PROBLEM ***\n");
  slap_FreeMatrix(S);
  return 0;
}
