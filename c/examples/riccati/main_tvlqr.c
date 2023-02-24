// Check README.md
// Sources: Lec. 7 Code on double integrator
//TODO: convert this to test Riccati

#include "slap/slap.h"
#include "riccati.h"

#define H 0.1
#define NSTATES 2
#define NINPUTS 1
#define NHORIZON 101

// array data to construct matrix, each column
double Q_data[NSTATES * NSTATES] = {1., 0., 0., 1.};
double R_data[NINPUTS * NINPUTS] = {0.1};
double q_data[NSTATES] = {0., 0.};
double r_data[NINPUTS] = {0.};

// A and B must be lists of matrices
double A_data[NSTATES * NSTATES * (NHORIZON - 1)];
double B_data[NSTATES * NINPUTS * (NHORIZON - 1)];

double xf_data[NSTATES] = {2., 0.};  //equilibrium point
double uf_data[NINPUTS] = {0.};
double x0_data[NSTATES] = {-1., 1.};

double Pp_data[NSTATES * (NSTATES + 1) * NHORIZON]; //stores P and p
double Kd_data[NINPUTS * (NSTATES + 1) * NHORIZON]; //stores K and d 

double x_data[NSTATES * NHORIZON];
double u_data[NINPUTS * (NHORIZON - 1)];
double y_data[NSTATES * NHORIZON];  // dual variables lambda

int main(void) {
  // Create matrix from array data
  Matrix Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  printf("\nQ = \n");
  slap_PrintMatrix(Q);
  Matrix R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  printf("\nR = ");
  slap_PrintMatrix(R);
  Matrix q = slap_MatrixFromArray(NSTATES, 1, q_data);
  Matrix r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  Matrix A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  printf("\nA = \n");
  slap_PrintMatrix(A);
  Matrix B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  printf("\nB = \n");
  slap_PrintMatrix(B);
  printf("End of problem\n");
  Matrix xf = slap_MatrixFromArray(NSTATES, 1, xf_data);
  Matrix uf = slap_MatrixFromArray(NINPUTS, 1, uf_data);
  Matrix x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);

  // Matrix of pointers
  Matrix Phist[NHORIZON];
  Matrix phist[NHORIZON];
  Matrix Khist[NHORIZON];
  Matrix dhist[NHORIZON];
  Matrix xhist[NHORIZON];
  Matrix uhist[NHORIZON];
  Matrix yhist[NHORIZON]; // dual variable: lambda

  // Pointer to the pre-allocated array
  double *Pp = Pp_data;
  double *Kd = Kd_data;
  double *xp = x_data;
  double *up = u_data;
  double *yp = y_data;
  for (int k = 0; k < NHORIZON; ++k) {
    // Pointer to each block, then next 
    Phist[k] = slap_MatrixFromArray(NSTATES, NSTATES, Pp);
    Pp += NSTATES * NSTATES;
    phist[k] = slap_MatrixFromArray(NSTATES, 1, Pp);
    Pp += NSTATES;
    Khist[k] = slap_MatrixFromArray(NINPUTS, NSTATES, Kd);
    Kd += NINPUTS * NSTATES;
    dhist[k] = slap_MatrixFromArray(NINPUTS, 1, Kd);
    Kd += NINPUTS;

    xhist[k] = slap_MatrixFromArray(NSTATES, 1, xp);
    xp += NSTATES;
    uhist[k] = slap_MatrixFromArray(NINPUTS, 1, up);
    up += NINPUTS;
    yhist[k] = slap_MatrixFromArray(NSTATES, 1, yp);
    yp += NSTATES;
  }
  // End of initializing memory and variables

  // Temporary matrix for underlying calculation 
  Matrix S = slap_NewMatrixZeros(NSTATES + NINPUTS, NSTATES + NINPUTS + 1);

  // Could work for delta_x and delta_u as well
  tiny_Riccati_LTV(NHORIZON - 1, A, B, Q, R, q, r, 
                   Khist, dhist, Phist, phist, S);
  // for (int k = 0; k < NHORIZON-1; k+=1) {
  //   printf("K[%d] = ", k);
  //   slap_PrintMatrix(Khist[k]);
  // }
  // tiny_LQR_LTI(NHORIZON - 1, A, B, Q, R, q, r, 
  //                  Khist, dhist, Phist, phist, S);
  tiny_RiccatiForwardPass_LTV(NHORIZON - 1, A, B, x0, xf, uf, 
                              Khist, dhist, Phist, phist, xhist, uhist, yhist);
  for (int k = 0; k < NHORIZON-1; k+=1) {
    printf("x[%d] = ", k);
    slap_PrintMatrix(slap_Transpose(xhist[k]));
    printf("u[%d] = ", k);
    slap_PrintMatrix(slap_Transpose(uhist[k]));
  }
  
  // printf("y[N] = ");
  // slap_PrintMatrix(slap_Transpose(yhist[NHORIZON-1]));
  slap_FreeMatrix(S);
  return 0;
}