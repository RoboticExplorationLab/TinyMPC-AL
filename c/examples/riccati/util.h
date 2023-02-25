#pragma once

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include "slap/slap.h"

int tiny_ReadData(const char* filename, double* des, const int size, bool verbose)
{
    FILE *input;
    int i;

    input = fopen(filename, "r");
    if (!input) {
        if (verbose == true)
            fprintf(stderr, "Cannot open %s: %s.\n", filename, strerror(errno));
        return EXIT_FAILURE;
    }

    for (i = 0; i < size; ++i) {
        if (fscanf(input, "%lf ", &(des[i])) != 1) {
            if (verbose == true)
                fprintf(stderr, "Invalid data in %s.\n", filename);
            fclose(input);
            return EXIT_FAILURE;
        }

        if (verbose == true) 
            printf("Read %lf from %s.\n", des[i], filename);
    }

    if (ferror(input)) {
        fclose(input);
        if (verbose == true)
            fprintf(stderr, "Error reading %s.\n", filename);
        return EXIT_FAILURE;
    }
    if (fclose(input)) {
        if (verbose == true)
            fprintf(stderr, "Error closing %s.\n", filename);
        return EXIT_FAILURE;
    }

    if (verbose == true)
        printf("All doubles read successfully.\n");

    return EXIT_SUCCESS;
}

// Matrix tiny_RK4(const Matrix (*dynamics)(const Matrix, const Matrix), 
//                 const Matrix x, const Matrix u, const float h)
// {
//     Matrix k1 = (*dynamics)(x, u);
//     slap_MatMulAdd(x, )
//     Matrix k2 = (*dynamics)(x + 0.5*h*k1, uk);
//     Matrix k3 = (*dynamics)(model, xk+0.5*h*k2, uk);
//     Matrix k4 = (*dynamics)(model, xk+h*k3, uk);
//     Matrix x = xk + h*(k1+2k2+2k3+k4)/6;
//     return x;
// }