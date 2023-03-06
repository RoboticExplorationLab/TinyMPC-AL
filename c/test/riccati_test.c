#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "augmented_lagrangian_lqr.h"
#include "simpletest/simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 3

void IneqInputsTest() {
  double u_max_data[2] = {1, 1};
  double u_min_data[2] = {-1, -1};
  double u_data[2] = {1.1, 0.9};

  tiny_ProblemData prob = kDefaultProblemData;
  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.ncstr_states = 2*NSTATES*NHORIZON;
  prob.ncstr_inputs = 2*NINPUTS*(NHORIZON-1);
  prob.ncstr_goal = NSTATES;
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
  Matrix u = slap_MatrixFromArray(NINPUTS, 1, u_data);
  Matrix ineq = tiny_IneqInputs(prob, u);
}

int main() {
  IneqInputsTest();
  PrintTestResult();
  return TestResult();
}
