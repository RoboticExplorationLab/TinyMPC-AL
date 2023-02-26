#pragma once

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>

#include "slap/slap.h"

//========================================
// Print matrix with its name (dummy)
//========================================
#define tiny_Print(mat) {printf("\n%s = \n", #mat); slap_PrintMatrix(mat);}

//========================================
// Read data from file
//========================================
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