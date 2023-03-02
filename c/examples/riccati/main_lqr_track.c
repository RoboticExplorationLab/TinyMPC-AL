// Check README.md
// Sources: Lec. 7 Code on double integrator
// Task: LTI system track a nominal trajectory in MPC style
// TODO: convert this to test Riccati

#include "riccati.h"
#include "slap/slap.h"
#include "util.h"

#define H 0.1
#define NSTATES 2
#define NINPUTS 1
#define NHORIZON 10
#define NSIM 101
// array data to construct matrix, each column
double Q_data[NSTATES * NSTATES] = {1., 0., 0., 1.};
double R_data[NINPUTS * NINPUTS] = {0.1};
double q_data[NSTATES] = {0., 0.};
double r_data[NINPUTS] = {0.};
double A_data[NSTATES * NSTATES] = {1, 0, H, 1};
double B_data[NSTATES * NINPUTS] = {0.5 * H * H, H};

double x0_data[NSTATES] = {2., 0.5};
double xn_data[NSTATES * NSIM];        // nominal states
double un_data[NINPUTS * (NSIM - 1)];  // nominal inputs

double Pp_data[NSTATES * (NSTATES + 1) * NHORIZON];  // stores P and p
double Kd_data[NINPUTS * (NSTATES + 1) * NHORIZON];  // stores K and d

double x_data[NSTATES * NSIM];
double u_data[NINPUTS * (NSIM - 1)];

int main(void) {
  // Read reference trajectory from files
  const char *file_xn = "../examples/riccati/data/xn_data.txt";
  const char *file_un = "../examples/riccati/data/un_data.txt";
  int size_xn = NSTATES * NSIM;
  int size_un = NINPUTS * (NSIM - 1);
  tiny_ReadData(file_xn, xn_data, size_xn, false);
  tiny_ReadData(file_un, un_data, size_un, false);

  // Create matrix from array data
  Matrix Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  Matrix R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  Matrix q = slap_MatrixFromArray(NSTATES, 1, q_data);
  Matrix r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  Matrix A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  Matrix B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  Matrix x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);

  // Matrix of pointers
  Matrix Phist[NHORIZON];
  Matrix phist[NHORIZON];
  Matrix Khist[NHORIZON];
  Matrix dhist[NHORIZON];
  Matrix xhist[NSIM];
  Matrix uhist[NSIM - 1];
  Matrix un[NSIM - 1];
  Matrix xn[NSIM];

  // Pointer to the pre-allocated array
  double *Pp = Pp_data;
  double *Kd = Kd_data;
  double *xp = x_data;
  double *up = u_data;
  double *unp = un_data;
  double *xnp = xn_data;

  for (int k = 0; k < NSIM; ++k) {
    // Pointer to each block, then next
    if (k < NHORIZON) {
      Phist[k] = slap_MatrixFromArray(NSTATES, NSTATES, Pp);
      Pp += NSTATES * NSTATES;
      phist[k] = slap_MatrixFromArray(NSTATES, 1, Pp);
      Pp += NSTATES;
      Khist[k] = slap_MatrixFromArray(NINPUTS, NSTATES, Kd);
      Kd += NINPUTS * NSTATES;
      dhist[k] = slap_MatrixFromArray(NINPUTS, 1, Kd);
      Kd += NINPUTS;
    }
    xhist[k] = slap_MatrixFromArray(NSTATES, 1, xp);
    xp += NSTATES;
    xn[k] = slap_MatrixFromArray(NSTATES, 1, xnp);
    xnp += NSTATES;

    if (k < NSIM - 1) {
      uhist[k] = slap_MatrixFromArray(NINPUTS, 1, up);
      up += NINPUTS;
      un[k] = slap_MatrixFromArray(NINPUTS, 1, unp);
      unp += NINPUTS;
    }
  }
  // End of initializing memory and variables

  // Temporary matrix for underlying calculation
  Matrix S = slap_NewMatrixZeros(NSTATES + NINPUTS, NSTATES + NINPUTS + 1);
  slap_MatrixCopy(xhist[0], x0);

  // MPC loop
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    tiny_LQR_LTI(NHORIZON - 1, A, B, Q, R, q, r, Khist, dhist, Phist, phist, S);

    // Control input: u = uf - d - K*(x - xf)
    slap_MatrixAddition(uhist[k], un[k], dhist[0], -1);   // u[k] = un[k] + d[k]
    slap_MatMulAdd(uhist[k], Khist[0], xhist[k], -1, 1);  // u[k] -= K[k] * x[k]
    slap_MatMulAdd(uhist[k], Khist[0], xn[k], 1, 1);  // u[k] += K[k] * xn[k]

    // Next state: x = A*x + B*u
    slap_MatMulAdd(xhist[k + 1], A, xhist[k], 1, 0);  // x[k+1] = A * x[k]
    slap_MatMulAdd(xhist[k + 1], B, uhist[k], 1, 1);  // x[k+1] += B * u[k]

    printf("ex[%d] = %.4f\n", k, slap_MatrixNormedDifference(xn[k], xhist[k]));
  }
  slap_FreeMatrix(S);
  return 0;
}