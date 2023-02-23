#include "slap_arduino.h"
#include <Wire.h>

double data_A[4];
double data_B[4];
double data_C[4];
double time = 0.;

void setup() {
/* serial to display data */
Serial.begin(115200);
while(!Serial) {}
  Matrix A = slap_MatrixFromArray(2, 2, data_A);
  Matrix B = slap_MatrixFromArray(2, 2, data_B);
  Matrix C = slap_MatrixFromArray(2, 2, data_C);
  slap_SetIdentity(A, 3);
  slap_SetConst(B, 2);
  Serial.println("A00");
  Serial.println(*slap_GetElement(A, 0, 0));
  Serial.println("B00");
  Serial.println(*slap_GetElement(B, 0, 0));
  // C = alpha * A * B + beta * C;
  double alpha = 2;
  double beta = 0;
  slap_MatMulAdd(C, A, B, alpha, beta);
  Serial.println("C00");
  Serial.println(*slap_GetElement(C, 0, 0));  
  Serial.println("B00");
  Serial.println(*slap_GetElement(B, 0, 0));    
  pinMode(13, OUTPUT);
}

void loop() {
  digitalWrite(13, HIGH);
  delay(500);
  digitalWrite(13, LOW);
  delay(500);
}
