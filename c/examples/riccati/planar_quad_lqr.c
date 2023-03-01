// Check README.md
// Sources: Lec. 9 Code on planar quadrotor
// New: Code generation for Jacobians and dynamics
// Task: Model is linearized about hovering. Run at random initial state and 
// go back to origin.
//TODO: convert this to test Riccati

#include <stdio.h>
#include "slap/slap.h"
#include "riccati.h"
#include "quadrotor.h"
#include "util.h"

#define H 0.1
#define NSTATES 6
#define NINPUTS 2
#define NHORIZON 40
#define NSIM (61 + NHORIZON)

int main(void) {
  printf("\n*** PROBLEM DEFINITION ***\n");
  struct tiny_Model_PlanarQuadrotor model = tiny_Model_Default;
  model.umin[0] = 0.2*model.m*model.g;
  model.umin[1] = model.umin[0];
  model.umax[0] = 0.6*model.m*model.g; model.umax[1] = model.umax[0];

  // array data to construct matrix, each column
  double Q_data[NSTATES * NSTATES] = {0};
  double R_data[NINPUTS * NINPUTS] = {0.1};
  double q_data[NSTATES] = {0};
  double r_data[NINPUTS] = {0};

  // A and B must be lists of matrices
  double A_data[NSTATES * NSTATES];
  double B_data[NSTATES * NINPUTS];

  double x0_data[NSTATES] = {2, 3, 0.2, 0.1, 0.1, 0.1};  // initial state
  double xref_data[NSTATES] = {0, 1, 0, 0, 0, 0};  // goal state
  double xhover_data[NSTATES] = {0};
  double uhover_data[NINPUTS] = {0.5*model.m*model.g, 0.5*model.m*model.g};

  double Pp_data[NSTATES * (NSTATES + 1) * NHORIZON]; //stores P and p
  double Kd_data[NINPUTS * (NSTATES + 1) * (NHORIZON - 1)]; //stores K and d 

  double x_data[NSTATES * NSIM];
  double u_data[NINPUTS * (NSIM - 1)];
  
  // Create matrix from array data
  Matrix Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(Q, 1);
  // printf("\nQ = \n"); slap_PrintMatrix(Q);
  tiny_Print(Q);  // test my macro
  Matrix R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(R, 0.01);
  printf("\nR = \n"); slap_PrintMatrix(R);
  Matrix q = slap_MatrixFromArray(NSTATES, 1, q_data);
  printf("\nq = \n"); slap_PrintMatrix(q);  
  Matrix r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  printf("\nr = \n"); slap_PrintMatrix(r);  
  Matrix A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  Matrix B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  
  Matrix xhover = slap_MatrixFromArray(NSTATES, 1, xhover_data);
  printf("\nxhover = \n"); slap_PrintMatrix(xhover);
  Matrix uhover = slap_MatrixFromArray(NINPUTS, 1, uhover_data);
  printf("\nuhover = \n"); slap_PrintMatrix(uhover);
  Matrix x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  printf("\nx0 = \n"); slap_PrintMatrix(x0);
  Matrix xref = slap_MatrixFromArray(NSTATES, 1, xref_data);
  tiny_Print(xref);
  // Matrix of pointers
  Matrix Phist[NHORIZON];
  Matrix phist[NHORIZON];
  Matrix Khist[NHORIZON - 1];
  Matrix dhist[NHORIZON - 1];
  Matrix xhist[NSIM];
  Matrix uhist[NSIM - 1];

  // Pointer to the pre-allocated array
  double *Pp = Pp_data;
  double *Kd = Kd_data;
  double *xp = x_data;
  double *up = u_data;

  for (int k = 0; k < NSIM; ++k) 
  {
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
    
    if (k < NSIM - 1)
    {
      uhist[k] = slap_MatrixFromArray(NINPUTS, 1, up);
      up += NINPUTS;
    }
  }
  // End of initializing memory and variables
  // Get hovering A and B
  tiny_GetJacobianA_Raw(A_data, xhover_data, uhover_data);
  tiny_GetJacobianB_Raw(B_data, xhover_data, uhover_data);
  printf("\nA = \n"); slap_PrintMatrix(A);
  printf("\nB = \n"); slap_PrintMatrix(B);
  // Temporary matrix for underlying calculation 
  Matrix S = slap_NewMatrixZeros(NSTATES + NINPUTS, NSTATES + NINPUTS + 1);
  slap_MatrixCopy(xhist[0], x0);  
  printf("\n*** START SOLVING ***\n");
  // MPC loop
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) 
  {
    tiny_LQR_LTI(NHORIZON - 1, A, B, Q, R, q, r, 
                     Khist, dhist, Phist, phist, S);

    // Control input: u = uf - d - K*(x - xf) 
    slap_MatrixAddition(uhist[k], uhover, dhist[0], -1);    // u[k] = un[k] - d[k]
    slap_MatMulAdd(uhist[k], Khist[0], xhist[k], -1, 1);   // u[k] -= K[k] * x[k]
    slap_MatMulAdd(uhist[k], Khist[0], xref, 1, 1);   // u[k] += K[k] * xn[k]

    // Next state: x = f(x, u)
    tiny_Clamps(uhist[k].data, model.umin, model.umax, NINPUTS);
    tiny_Dynamics_RK4(xhist[k+1], xhist[k], uhist[k]);
    printf("ex[%d] = %.4f\n", k, slap_MatrixNormedDifference(xref, xhist[k]));
    // printf("x[%d] = ", k);
    // slap_PrintMatrix(slap_Transpose(xhist[k]));
    // printf("u[%d] = ", k);
    // slap_PrintMatrix(slap_Transpose(uhist[k]));
  }
  slap_FreeMatrix(S);
  return 0;
}
