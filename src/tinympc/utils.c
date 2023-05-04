#include "utils.h"

void PrintSolveInfo(tiny_Workspace* work) {
  tiny_Info* info = work->info;
  printf("Solve info: \n");
  printf(" Status: %d\n", info->status_val);
  printf(" Iter Riccati: %d, Iter AL: %d\n", info->iter_riccati, info->iter_al);   
  printf(" Primal obj: %f, AL obj: %f\n", info->obj_pri, info->obj_al);
  printf(" Primal res: %f, dual res: %f\n", info->pri_res, info->dua_res);
}


// //========================================
// // Read data from file
// //========================================
// int tiny_ReadData(const char* filename, sfloat* des, const int size,
//                   bool verbose) {
//   FILE* input;
//   int i;

//   input = fopen(filename, "r");
//   if (!input) {
//     if (verbose == true)
//       fprintf(stderr, "Cannot open %s: %s.\n", filename, strerror(errno));
//     return EXIT_FAILURE;
//   }

//   for (i = 0; i < size; ++i) {
//     if (fscanf(input, "%lf ", &(des[i])) != 1) {
//       if (verbose == true) fprintf(stderr, "Invalid data in %s.\n", filename);
//       fclose(input);
//       return EXIT_FAILURE;
//     }

//     if (verbose == true) printf("Read %lf from %s.\n", des[i], filename);
//   }

//   if (ferror(input)) {
//     fclose(input);
//     if (verbose == true) fprintf(stderr, "Error reading %s.\n", filename);
//     return EXIT_FAILURE;
//   }
//   if (fclose(input)) {
//     if (verbose == true) fprintf(stderr, "Error closing %s.\n", filename);
//     return EXIT_FAILURE;
//   }

//   if (verbose == true) printf("All sfloats read successfully.\n");

//   return EXIT_SUCCESS;
// }

// //========================================
// // Read data from file and copy the last knot point into
// // remaining space of the array. Useful for extend horizon at the end.
// //========================================
// int tiny_ReadData_Extend(const char* filename, sfloat* des, const int stride,
//                          const int size, bool verbose) {
//   FILE* input;
//   int i;
//   int k = 0;
//   input = fopen(filename, "r");
//   if (!input) {
//     if (verbose == true)
//       fprintf(stderr, "Cannot open %s: %s.\n", filename, strerror(errno));
//     return EXIT_FAILURE;
//   }

//   for (i = 0; i < size; ++i) {
//     if (fscanf(input, "%lf ", &(des[i])) != 1) {
//       if (verbose == true) fprintf(stderr, "Invalid data in %s.\n", filename);
//       fclose(input);
//       break;
//     }

//     if (verbose == true) printf("Read %lf from %s.\n", des[i], filename);

//     k += 1;
//   }

//   if (verbose == true)
//     printf("All sfloats read successfully and now extend.\n");

//   int remain_cnt = (size - k) / stride;  // # of remaining chunks
//   for (i = 0; i < remain_cnt; i += 1) {
//     for (int j = 0; j < stride; j += 1) {
//       des[k + j + i * stride] = des[k + j - stride];  // copy
//     }
//   }

//   return EXIT_SUCCESS;
// }

// //========================================
// // Read data from file and copy the goal state into
// // remaining space of the array. Useful for extend horizon at the end.
// //========================================
// int tiny_ReadData_ExtendGoal(const char* filename, sfloat* des,
//                              const sfloat* xf, const int stride, const int size,
//                              bool verbose) {
//   FILE* input;
//   int i;
//   int k = 0;
//   input = fopen(filename, "r");
//   if (!input) {
//     if (verbose == true)
//       fprintf(stderr, "Cannot open %s: %s.\n", filename, strerror(errno));
//     return EXIT_FAILURE;
//   }

//   for (i = 0; i < size; ++i) {
//     if (fscanf(input, "%lf ", &(des[i])) != 1) {
//       if (verbose == true) fprintf(stderr, "Invalid data in %s.\n", filename);
//       fclose(input);
//       break;
//     }

//     if (verbose == true) printf("Read %lf from %s.\n", des[i], filename);

//     k += 1;
//   }

//   if (verbose == true)
//     printf("All sfloats read successfully and now extend.\n");

//   int remain_cnt = (size - k) / stride;  // # of remaining chunks
//   for (i = 0; i < remain_cnt; i += 1) {
//     for (int j = 0; j < stride; j += 1) {
//       des[k + j + i * stride] = xf[j];  // copy
//     }
//   }

//   return EXIT_SUCCESS;
// }

//========================================
// Clamp the inputs to within min max value
//========================================
void tiny_Clamps(sfloat* arr, const sfloat* min, const sfloat* max,
                 const int N) {
  for (int k = 0; k < N; ++k) {
    arr[k] = (arr[k] > max[k]) ? max[k] : ((arr[k] < min[k]) ? min[k] : arr[k]);
  }
}

void tiny_Clamp(sfloat* arr, const sfloat min, const sfloat max, const int N) {
  for (int k = 0; k < N; ++k) {
    arr[k] = (arr[k] > max) ? max : ((arr[k] < min) ? min : arr[k]);
  }
}

// Clamp all data for matrix or vector
void tiny_ClampMatrix(Matrix* mat, const Matrix min, const Matrix max) {
  tiny_Clamps(mat->data, min.data, max.data, (mat->rows) * (mat->cols));
}

void tiny_ShiftFill(Matrix* mats, const int length) {
  for (int k = 0; k < length - 1; ++k) {
    slap_Copy(mats[k], mats[k + 1]);
  }
  slap_Copy(mats[length - 1], mats[length - 2]);
}

void tiny_ShiftFillWith(Matrix* mats, const sfloat* x, const int length) {
  for (int k = 0; k < length - 1; ++k) {
    slap_Copy(mats[k], mats[k + 1]);
  }
  slap_CopyFromArray(mats[length - 1], x);
}