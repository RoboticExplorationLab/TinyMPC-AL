#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "augmented_lagrangian_lqr.h"
#include "simpletest/simpletest.h"
// #include "test_utils.h"

void MyTest() {
    
}

int main() {
  TEST(2 == 1 + 1);
  TEST(2 == 2);
  MyTest();
  PrintTestResult();
  return TestResult();
}
