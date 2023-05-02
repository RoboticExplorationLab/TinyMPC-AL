#include "types.h"

sfloat tiny_RiccatiConvergence(const tiny_ProblemData prob);

void tiny_IneqInputs(Matrix* ineq, const tiny_ProblemData prob, const Matrix u);

// void tiny_IneqInputsOffset(Matrix* cu, const tiny_ProblemData prob);

// void tiny_IneqInputsJacobian(Matrix* ineq_jac, const tiny_ProblemData prob);

void tiny_IneqStates(Matrix* cx, const tiny_ProblemData prob,
                     const Matrix x);

// void tiny_IneqStatesOffset(Matrix* cx, const tiny_ProblemData prob);

// void tiny_IneqStatesJacobian(Matrix* ineq_jac, const tiny_ProblemData prob);

void tiny_ActiveIneqMask(Matrix* mask, const Matrix input_dual,
                         const Matrix ineq);

void tiny_ClampIneqDuals(Matrix* dual, const Matrix new_dual);