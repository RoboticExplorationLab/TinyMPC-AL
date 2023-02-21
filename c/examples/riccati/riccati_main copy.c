//
// Created by Brian Edward Jackson on 1/31/23.
//

#include "slap/slap.h"
#include "riccati.h"

#define NSTATES 3
#define NINPUTS 2
#define NHORIZON 5

// array data to construct matrix
double Q_data[NSTATES * NSTATES] = {1, 0, 0,  
                                    0, 0.5, 0,  
                                    0, 0, 0.4};
double R_data[NINPUTS * NINPUTS] = {0.1, 0,  0, 0.2};
double q_data[NSTATES] = {-0.1, -0.2, -0.3};
double r_data[NINPUTS] = {0.1, 0.2};

double A_data[NSTATES * NSTATES] = {1, 0, 0,  
                                    0, 1, 0,  
                                    0, 0, 1};
double B_data[NSTATES * NINPUTS] = {.1, 0,  
                                    0, 0,  
                                    .1, 0};
double f_data[NSTATES] = {-0.1, 0.2, 0.3};
double x0_data[NSTATES] = {1,2,3};

double Pp_data[NSTATES * (NSTATES + 1) * NHORIZON];
double Kd_data[NINPUTS * (NSTATES + 1) * NHORIZON];

double x_data[NSTATES * NHORIZON];
double u_data[NINPUTS * (NHORIZON - 1)];
double y_data[NSTATES * NHORIZON];  // dual variables

int main(void) {
  // Create matrix from array data
  Matrix Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  Matrix R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  Matrix q = slap_MatrixFromArray(NSTATES, 1, q_data);
  Matrix r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  Matrix A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  Matrix B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  Matrix f = slap_MatrixFromArray(NSTATES, 1, f_data);
  Matrix x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);

  // Matrix of pointers
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];
  Matrix K[NHORIZON];
  Matrix d[NHORIZON];
  Matrix x[NHORIZON];
  Matrix u[NHORIZON];
  Matrix y[NHORIZON];

  // Pointer to the pre-allocated array
  double *Pp = Pp_data;
  double *Kd = Kd_data;
  double *xp = x_data;
  double *up = u_data;
  double *yp = y_data;
  for (int k = 0; k < NHORIZON; ++k) {
    // Pointer to each block, then next
    P[k] = slap_MatrixFromArray(NSTATES, NSTATES, Pp);
    Pp += NSTATES * NSTATES;
    p[k] = slap_MatrixFromArray(NSTATES, 1, Pp);
    Pp += NSTATES;
    K[k] = slap_MatrixFromArray(NINPUTS, NSTATES, Kd);
    Kd += NINPUTS * NSTATES;
    d[k] = slap_MatrixFromArray(NINPUTS, 1, Kd);
    Kd += NINPUTS;

    x[k] = slap_MatrixFromArray(NSTATES, 1, xp);
    xp += NSTATES;
    u[k] = slap_MatrixFromArray(NINPUTS, 1, up);
    up += NINPUTS;
    y[k] = slap_MatrixFromArray(NSTATES, 1, yp);
    yp += NSTATES;
  }
  // End of initializing memory and variables

  // Temporary matrix for underlying calculation 
  Matrix S = slap_NewMatrixZeros(NSTATES + NINPUTS, NSTATES + NINPUTS + 1);
  slap_Riccati_LTI(NHORIZON - 1, A, B, f, Q, R, q, r, K, d, P, p, S);

  slap_RiccatiForwardPass_LTI(NHORIZON - 1, A, B, f, x0, K, d, P, p, x, u, y);
  printf("x[N] = ");
  slap_PrintMatrix(slap_Transpose(x[NHORIZON-1]));
  printf("u[N-1] = ");
  slap_PrintMatrix(slap_Transpose(u[NHORIZON-2]));
  printf("y[N] = ");
  slap_PrintMatrix(slap_Transpose(y[NHORIZON-1]));
  slap_FreeMatrix(S);
  return 0;
}
