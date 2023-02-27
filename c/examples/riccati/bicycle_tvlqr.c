// Check README.md
// Second-order bicycle model tracking TVLQR
// Code generation for Jacobians and dynamics
//TODO: convert this to test Riccati

#include <stdio.h>
#include "util.h"
#include "slap/slap.h"
#include "riccati.h"
#include "bicycle.h"

#define H 0.1
#define NSTATES 5
#define NINPUTS 2
#define NHORIZON 30
#define NSIM (101 + 1*NHORIZON)  // extend the end to complete MPC

int main(void) {
  printf("\n*** PROBLEM DEFINITION ***\n");
  struct tiny_Model_Bicycle model;
  model.x45max[0] = 2;      // velocity max
  model.x45max[1] = 0.5;    // steering rate max
  model.umax[0] = 4;        // accel max
  model.umax[1] = 0.7;      // steering angle max
  model.x45min[0] = -model.x45max[0]; model.x45min[1] = -model.x45max[1];
  model.umin[0] = -model.umax[0];model.umin[1] = -model.umax[1];

  // array data to construct matrix, each column
  double Q_data[NSTATES * NSTATES] = {0};
  double R_data[NINPUTS * NINPUTS] = {0};
  double q_data[NSTATES] = {0};
  double r_data[NINPUTS] = {0};

  // A and B must be lists of matrices
  double A_data[NSTATES * NSTATES] = {0};
  double B_data[NSTATES * NINPUTS] = {0};

  double x0_data[NSTATES] = {2, 0, 0, 0, 0};
  double xf_data[NSTATES] = {20, 20, M_PI/4, 0, 0};
  double uf_data[NINPUTS] = {0};

  double Pp_data[NSTATES * (NSTATES + 1) * NHORIZON] = {0}; //stores P and p
  double Kd_data[NINPUTS * (NSTATES + 1) * (NHORIZON - 1)] = {0}; //stores K and d 

  double x_data[NSTATES * NSIM] = {0};
  double u_data[NINPUTS * (NSIM - 1)] = {0};
  double xref_data[NSTATES * NSIM] = {0};
  double uref_data[NINPUTS * (NSIM - 1)] = {0}; 

  // Read reference trajectory from files
  const char *file_xref = "../examples/riccati/data/xref_data.txt";
  const char *file_uref = "../examples/riccati/data/uref_data.txt";
  tiny_ReadData(file_xref, xref_data, NSTATES * NSIM, false);
  tiny_ReadData(file_uref, uref_data, NINPUTS * (NSIM - 1), false);  
  tiny_ReadData_ExtendGoal(file_xref, xref_data, xf_data, NSTATES, NSTATES * NSIM, false);
  tiny_ReadData_ExtendGoal(file_uref, uref_data, uf_data, NINPUTS, NINPUTS * (NSIM-1), false);
  // tiny_ReadData_Extend(file_xref, xref_data, NSTATES, NSTATES * NSIM, false);
  // tiny_ReadData_Extend(file_uref, uref_data, NINPUTS, NINPUTS * (NSIM-1), false);  
  // Create matrix from array data
  Matrix Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(Q, 1);
  // printf("\nQ = \n"); slap_PrintMatrix(Q);
  tiny_Print(Q);
  Matrix R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(R, 0.01);
  printf("\nR = \n"); slap_PrintMatrix(R);
  Matrix q = slap_MatrixFromArray(NSTATES, 1, q_data);
  // printf("\nq = \n"); slap_PrintMatrix(q);  
  Matrix r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  // printf("\nr = \n"); slap_PrintMatrix(r);  
  Matrix A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  Matrix B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);

  Matrix x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  printf("\nx0 = \n"); slap_PrintMatrix(x0);
  Matrix xf = slap_MatrixFromArray(NSTATES, 1, xf_data);
  tiny_Print(xf);  

  // Matrix of pointers
  Matrix Phist[NHORIZON];
  Matrix phist[NHORIZON];
  Matrix Khist[NHORIZON - 1];
  Matrix dhist[NHORIZON - 1];
  Matrix xhist[NSIM];
  Matrix uhist[NSIM - 1];
  Matrix uref[NSIM - 1];
  Matrix xref[NSIM];
  // Pointer to the pre-allocated array
  double *Pp = Pp_data;
  double *Kd = Kd_data;
  double *xp = x_data;
  double *up = u_data;
  double *xrefp = xref_data;
  double *urefp = uref_data;

  for (int k = 0; k < NSIM; ++k) {
    // Pointer to each block, then next 
    if (k < NHORIZON)
    {
      Phist[k] = slap_MatrixFromArray(NSTATES, NSTATES, Pp);
      Pp += NSTATES * NSTATES;
      phist[k] = slap_MatrixFromArray(NSTATES, 1, Pp);
      Pp += NSTATES;
      if (k < NHORIZON - 1)
      {
        Khist[k] = slap_MatrixFromArray(NINPUTS, NSTATES, Kd);
        Kd += NINPUTS * NSTATES;
        dhist[k] = slap_MatrixFromArray(NINPUTS, 1, Kd);
        Kd += NINPUTS;
      }
    }
    xhist[k] = slap_MatrixFromArray(NSTATES, 1, xp);
    xp += NSTATES;
    xref[k] = slap_MatrixFromArray(NSTATES, 1, xrefp);
    xrefp += NSTATES;
    
    if (k < NSIM - 1)
    {
      uhist[k] = slap_MatrixFromArray(NINPUTS, 1, up);
      up += NINPUTS;
      uref[k] = slap_MatrixFromArray(NINPUTS, 1, urefp);
      urefp += NINPUTS;
    }
  }
  // End of initializing memory and variables

  // Temporary matrix for underlying calculation 
  Matrix S = slap_NewMatrixZeros(NSTATES + NINPUTS, NSTATES + NINPUTS + 1);
  slap_MatrixCopy(xhist[0], x0);  
  printf("\n*** START SOLVING ***\n");
  // MPC loop
  for (int k = 0; k < NSIM - NHORIZON; ++k) {
    tiny_LQR_LTVf(NHORIZON - 1, A, B, tiny_GetJacobians, Q, R, q, r, 
                     Khist, dhist, Phist, phist, &xref[k], &uref[k], S);
    // tiny_Print(Khist[0]);

    // Control input: u = uf - d - K*(x - xf) 
    slap_MatrixAddition(uhist[k], uref[k], dhist[0], -1);    // u[k] = un[k] - d[k]
    // slap_PrintMatrix(slap_Transpose(uhist[k]));
    slap_MatMulAdd(uhist[k], Khist[0], xhist[k], -1, 1);   // u[k] -= K[k] * x[k]
    // slap_PrintMatrix(slap_Transpose(uhist[k]));
    slap_MatMulAdd(uhist[k], Khist[0], xref[k], 1, 1);   // u[k] += K[k] * xn[k]
    // printf("u[%d] = ", k);
    // slap_PrintMatrix(slap_Transpose(uhist[k]));

    // Next state: x = f(x, u)
    tiny_Clamps(uhist[k].data, model.umin, model.umax, NINPUTS);
    tiny_Dynamics_RK4_Raw(xhist[k+1].data, xhist[k].data, uhist[k].data);

    // Tracking errors (show different metrics)
    printf("ex[%d] = %.4f\n", k, slap_MatrixNormedDifference(xref[k], xhist[k]));
    // printf("x[%d] = ", k);
    // slap_PrintMatrix(slap_Transpose(xhist[k]));

  }
  printf("\n*** END OF PROBLEM ***\n");
  slap_FreeMatrix(S);
  return 0;
}
