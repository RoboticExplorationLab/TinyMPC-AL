#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "cjson/cJSON.h" // or #include <cJSON/cJSON.h> since it's now an external library linked by cmake


//========================================
// Read data from file
//========================================
int read_data(const char* filename, double* dst, const int dst_size, bool verbose) {
  int res = 0;
  // Open file
  FILE* file = fopen(filename, "r");
  // Get file size
  fseek(file, 0L, SEEK_END);
  long int file_len = ftell(file);
  char file_buf[file_len];
  // Seek back to beginning of file
  rewind(file);
  // Write file char data into file_buf
  char c;
  for (long int i=0; i<file_len; i++) {
    c = fgetc(file);
    if (c == EOF) {
      res = 1;
    }
    file_buf[i] = c;
  }
  // Check file was read successfully
  if (res == 1) {
    if (verbose) {
      fprintf(stderr, "Failed to read file in read_data\n");
    }
  }
  // Add null terminator to file_buf
  file_buf[file_len] = '\0';

  // Print file contents (for debugging)
  // printf("file_buf: %s\n", file_buf);

  cJSON* root = cJSON_Parse(file_buf);
  cJSON* temp = root;

  // Populate destination with data from JSON
  // Point arr_ptr to first array in set of arrays
  //               root->outer->first arr
  cJSON* arr_ptr;
  cJSON* subarr_ptr;
  arr_ptr = temp->child->child;
  long int i = 0;
  while (arr_ptr) {
    // Get every value in array or subarray
    if (arr_ptr->child) {
      subarr_ptr = arr_ptr->child;
      while (subarr_ptr) {
        dst[i] = subarr_ptr->valuedouble;
        subarr_ptr = subarr_ptr->next;
        i++;
      }
    } else {
      char* jsonstr = cJSON_Print(arr_ptr);
      printf("jsonstr: %s\n", jsonstr);
      printf("valuedouble: %f\n", arr_ptr->valuedouble);
      dst[i] = arr_ptr->valuedouble;
      i++;
    }

    
    // Point arr_ptr to next array in set of arrays
    arr_ptr = arr_ptr->next;
  }

  if (i != dst_size) {
    fprintf(stderr, "Data available in json not equal to destination array size.\n");
    res = 1;
    return res;
  }

  // Delete cJSON pointers
  cJSON_Delete(root);
  cJSON_Delete(arr_ptr);
  cJSON_Delete(subarr_ptr);

  return res;
}

int main(void) {
  
  puts("we doin json stuff\n");

  const int NSTATES = 2;
  const int NINPUTS = 1;
  const int NSIM = 101;

  double xn_data[NSTATES * NSIM];        // nominal states
  double un_data[NINPUTS * (NSIM - 1)];  // nominal inputs

  int res = read_data("../examples/json/data/xn_data.json", xn_data, NSTATES * NSIM, true);
  if (res == EXIT_FAILURE) {
    fprintf(stderr, "Error when reading xn_data\n");
    return EXIT_FAILURE;
  }
  res = read_data("../examples/json/data/un_data.json", un_data, NINPUTS * (NSIM-1), true);
  if (res == EXIT_FAILURE) {
    fprintf(stderr, "Error when reading un_data\n");
    return EXIT_FAILURE;
  }
  
  printf("xn_data:\n");
  for (long int i=0; i<(int)(sizeof(xn_data)/sizeof(xn_data[0])); i++) {
    printf("%f\n", xn_data[i]);
  }
  printf("\n----------------------------------------------\n");
  printf("\nun_data:\n");
  for (long int i=0; i<(int)(sizeof(un_data)/sizeof(un_data[0])); i++) {
    printf("%f\n", un_data[i]);
  }

  return EXIT_SUCCESS;
}