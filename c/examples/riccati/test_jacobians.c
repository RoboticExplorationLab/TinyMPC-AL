// Check README.md
// Test Jacobians from generated functions

#include "slap/slap.h"
#include "riccati.h"
#include <stdio.h>
#include "util.h"

#define H 0.1
#define NSTATES 6
#define NINPUTS 2
#define NHORIZON 101

// A and B must be lists of matrices
double A_data[NSTATES * NSTATES];
double B_data[NSTATES * NINPUTS];

double x0_data[NSTATES] = {1, 2, 3, 4, 5, 6};
double u0_data[NINPUTS] = {7, 8};

double x_data[NSTATES * NHORIZON];
double u_data[NINPUTS * (NHORIZON - 1)];

int main(void) {
  Matrix A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  Matrix B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  printf("End of problem\n");
  // Matrix xf = slap_MatrixFromArray(NSTATES, 1, xf_data);
  // Matrix uf = slap_MatrixFromArray(NINPUTS, 1, uf_data);
  Matrix x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  Matrix u0 = slap_MatrixFromArray(NINPUTS, 1, u0_data);

  tiny_GetJacobianA(A_data, x0_data, u0_data);
  tiny_GetJacobianB(B_data, x0_data, u0_data);
  // Compare with results from julia code
  slap_PrintMatrix(A); 
  slap_PrintMatrix(B); 
  return 0;
}