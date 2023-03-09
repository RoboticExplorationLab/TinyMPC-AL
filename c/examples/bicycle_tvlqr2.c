// Check README.md
// Second-order bicycle model tracking TVLQR, test forward
// New: Code generation for Jacobians and dynamics
// Task: LTV system to track a reference trajectory
// TODO: convert this to test Riccati

#include <stdio.h>

#include "bicycle.h"
#include "riccati.h"
#include "slap/slap.h"
#include "util.h"

#define H 0.1
#define NSTATES 5
#define NINPUTS 2
#define NHORIZON 101

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

  double x_data[NSTATES * NHORIZON] = {0};
  double u_data[NINPUTS * (NHORIZON - 1)] = {0};
  double xref_data[NSTATES * NHORIZON] = {0};
  double uref_data[NINPUTS * (NHORIZON - 1)] = {0};

  // Read reference trajectory from files
  const char *file_xref = "../examples/riccati/data/xref_data.txt";
  const char *file_uref = "../examples/riccati/data/uref_data.txt";
  tiny_ReadData(file_xref, xref_data, NSTATES * NHORIZON, false);
  tiny_ReadData(file_uref, uref_data, NINPUTS * (NHORIZON - 1), false);
  // Create matrix from array data
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
  Matrix Phist[NHORIZON];
  Matrix phist[NHORIZON];
  Matrix Khist[NHORIZON - 1];
  Matrix dhist[NHORIZON - 1];
  Matrix xhist[NHORIZON];
  Matrix uhist[NHORIZON - 1];
  Matrix uref[NHORIZON - 1];
  Matrix xref[NHORIZON];
  // Pointer to the pre-allocated array
  double *Pp = Pp_data;
  double *Kd = Kd_data;
  double *xp = x_data;
  double *up = u_data;
  double *xrefp = xref_data;
  double *urefp = uref_data;

  for (int k = 0; k < NHORIZON; ++k) {
    // Pointer to each block, then next
    Phist[k] = slap_MatrixFromArray(NSTATES, NSTATES, Pp);
    Pp += NSTATES * NSTATES;
    phist[k] = slap_MatrixFromArray(NSTATES, 1, Pp);
    Pp += NSTATES;
    xhist[k] = slap_MatrixFromArray(NSTATES, 1, xp);
    xp += NSTATES;
    xref[k] = slap_MatrixFromArray(NSTATES, 1, xrefp);
    xrefp += NSTATES;

    if (k < NHORIZON - 1) {
      uhist[k] = slap_MatrixFromArray(NINPUTS, 1, up);
      up += NINPUTS;
      uref[k] = slap_MatrixFromArray(NINPUTS, 1, urefp);
      urefp += NINPUTS;
      Khist[k] = slap_MatrixFromArray(NINPUTS, NSTATES, Kd);
      Kd += NINPUTS * NSTATES;
      dhist[k] = slap_MatrixFromArray(NINPUTS, 1, Kd);
      Kd += NINPUTS;
    }
  }
  // End of initializing memory and variables

  // Temporary matrix for underlying calculation
  Matrix S = slap_NewMatrixZeros(NSTATES + NINPUTS, NSTATES + NINPUTS + 1);
  slap_MatrixCopy(xhist[0], x0);
  printf("\n*** START SOLVING ***\n");
  int test = 100 - 61;
  tiny_Riccati_LTVf(NHORIZON - test, A, B, tiny_GetJacobians, Q, R, q, r, Khist,
                    dhist, Phist, phist, xref, uref, S);
  tiny_RiccatiForwardPass_LTVf(NHORIZON - test, A, B, tiny_GetJacobians, x0,
                               xref, uref, Khist, dhist, Phist, phist, xhist,
                               uhist, NULL);
  for (int k = 0; k < NHORIZON - test; k += 1) {
    // Tracking errors (show different metrics)
    printf("ex[%d] = %.4f\n", k,
           slap_MatrixNormedDifference(xref[k], xhist[k]));
    // printf("x[%d] = ", k);
    // slap_PrintMatrix(slap_Transpose(uhist[k]));
    // tiny_Print(Khist[k]);
  }

  printf("\n*** END OF PROBLEM ***\n");
  slap_FreeMatrix(S);
  return 0;
}
