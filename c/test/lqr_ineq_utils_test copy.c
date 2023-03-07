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

double A_data[NSTATES*NSTATES] = {1,0,0,0, 0,1,0,0, 0.1,0,0.1,0, 0,0.1,0,0.1};
double B_data[NSTATES*NINPUTS] = {0.005,0,0.1,0, 0,0.005,0,0.1};
double f_data[NSTATES] = {0,0,0,0};
double x0_data[NSTATES] = {5,7,2,-1.4};

void ForwardPassTest() {
  tiny_LinearDiscreteModel model = kDefaultLinearDiscreteModel;
  model.dt = 0.1;
  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  model.B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  model.f = slap_MatrixFromArray(NSTATES, 1, f_data);
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);

}



int main() {
  ForwardPassTest();
  return TestResult();
}
