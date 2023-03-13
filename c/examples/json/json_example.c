#include <stdio.h>
// #include <stdbool.h>
// #include <errno.h>

// #include "cjson/cJSON.h" // or #include <cJSON/cJSON.h> since it's now an external library linked by cmake

// //========================================
// // Read data from file
// //========================================
// int read_data(const char* filename, double* des, const int size,
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

//   return EXIT_SUCCESS;
// }

int main(void) {
    puts("we doin json stuff\n");
}