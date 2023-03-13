#include "qpal.h"

enum slap_ErrorCode get_augmented_lagrangian_cost(const Matrix A, const Matrix b,
                            const Matrix Q, const Matrix R,
                            Matrix *x, Matrix *u, int *L, Matrix *L_grad, Matrix *L_hess,
                            const int NSTATES, const int NINPUTS, const int NHORIZON,
                            Matrix xf, Matrix q, Matrix r, Matrix *Khist,
                            Matrix *dhist, Matrix *Phist, Matrix *phist, Matrix *xhist,
                            Matrix *uhist, const Matrix yhist, Matrix x0, const int rho){
    
    for (int k = 0; k < NHORIZON - 1; k+=1){
        // (1/2)*(x_k - x_ref)' * Q * (x_k - x_ref)
        Matrix Q_matr;
        slap_MatMulAtB(Q_matr, x[k], Q);
        slap_MatMulAdd(Q_matr, Q_matr, x[k], 1, 0);

        // (1/2)*(u_k - u_ref)' * R * (u_k - u_ref)
        Matrix R_matr;
        slap_MatMulAtB(R_matr, u[k], R);
        slap_MatMulAdd(R_matr, R_matr, u[k], 1, 0);

        *L += 0.5 * (*Q_matr.data + *R_matr.data); // this is an integer
    }
    Matrix Qh_matr;
    slap_MatMulAtB(Qh_matr, xf, Q);
    slap_MatMulAdd(Qh_matr, Qh_matr, xf, 1, 0);
    *L += 0.5 * *Qh_matr.data; // this is an integer
    
}