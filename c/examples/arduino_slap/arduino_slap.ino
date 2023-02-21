#include "slap_arduino.h"

double data_A[12];
double data_B[6];
double data_C[8];

void setup() {
  Matrix A = slap_MatrixFromArray(4, 3, data_A);
  Matrix B = slap_MatrixFromArray(3, 2, data_B);
  Matrix C = slap_MatrixFromArray(4, 2, data_B);
  slap_SetConst(A, 5);
  slap_SetIdentity(B, 1.0);

  // C = alpha * A * B + beta * C;
  double alpha = 1;
  double beta = 0;
  slap_MatMulAdd(C, A, B, alpha, beta);

  pinMode(13, OUTPUT);
}

void loop() {
  digitalWrite(13, HIGH);
  delay(400);
  digitalWrite(13, LOW);
  delay(400);
}

