#include "tiny_cost.h"

void tiny_AddStageCost(double* cost, const tiny_ProblemData prob,
                       const Matrix x, const Matrix u, const int k) {
  double dx_data[prob.nstates];
  Matrix dx = slap_MatrixFromArray(prob.nstates, 1, dx_data);
  slap_MatrixAddition(dx, x, prob.X_ref[k], -1);
  *cost += 0.5 * slap_QuadraticForm(dx, prob.Q, dx);
  Matrix du = slap_MatrixFromArray(prob.ninputs, 1, dx_data);
  slap_MatrixAddition(du, u, prob.U_ref[k], -1);
  *cost += 0.5 * slap_QuadraticForm(du, prob.R, du);
}

void tiny_AddTerminalCost(double* cost, const tiny_ProblemData prob,
                          const Matrix x) {
  double dx_data[prob.nstates];
  Matrix dx = slap_MatrixFromArray(prob.nstates, 1, dx_data);
  slap_MatrixAddition(dx, x, prob.X_ref[prob.nhorizon - 1], -1);
  *cost += 0.5 * slap_QuadraticForm(dx, prob.Qf, dx);
}

void tiny_ExpandStageCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                          Matrix* hes_el_uu, Matrix* grad_el_u,
                          const tiny_ProblemData prob, const int k) {
  slap_MatrixCopy(*hes_el_xx, prob.Q);
  slap_MatMulAdd(*grad_el_x, prob.Q, prob.X_ref[k], -1, 0);
  slap_MatrixCopy(*hes_el_uu, prob.R);
  slap_MatMulAdd(*grad_el_u, prob.R, prob.U_ref[k], -1, 0);
}

void tiny_ExpandTerminalCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                             const tiny_ProblemData prob) {
  slap_MatrixCopy(*hes_el_xx, prob.Qf);
  slap_MatMulAdd(*grad_el_x, prob.Qf, prob.X_ref[prob.nhorizon-1], -1, 0);
}