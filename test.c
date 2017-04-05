#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/*
this is just a test file for me to test and mess around with code
*/

// mat0 is product, mat1, mat2 are things to be multiplied
void mat_mul(long double prod[3][3], long double mat1[3][3], long double mat2[3][3]) {
  int i, j, k;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      prod[i][j] = 0.0;
      for (k = 0; k < 3; k++) {
        prod[i][j] += mat1[i][k] * mat2[k][j];
      }
    }
  }

  return; 
}

// mat0 stores the product, mat2 is the thing to be exponentiated
void mat_pow(long double mat0[3][3], long double mat2[3][3], int power) {
  long double mat1[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        mat1[i][j] = mat2[i][j];
      }
  }

  // if power == 0, set mat0 to identity matrix
  if (power == 0) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
          if (i == j) {
            mat0[i][j] = 1.0;
          }
          else {
            mat0[i][j] = 0.0;
          }
      }
    }
  }
  
  // if power == 1, set mat0 to mat2
  else if (power == 1) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          mat0[i][j] = mat2[i][j];
        }
      }
  }

  for (int p = 1; p < power; p++) {
    // square mat2 and let mat0 equal mat2**2?
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          mat0[i][j] = 0.0;
          for (int k = 0; k < 3; k++) {
              mat0[i][j] += mat1[i][k]*mat2[k][j];
            }
        }
      }
    // let mat1 = mat2**2
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          mat1[i][j] = mat0[i][j];
        }
      }
  }

  return; 
}

void print_m(long double m[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%Lf ", m[i][j]);
        }
        printf("\n");
    }
}

void r_mat_pow (long double prod[3][3], long double mat[3][3], int power) {
  int i, j, k;

  // set to identity matrix
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      if (i == j) {
        prod[i][j] = 1.0;
      }
      else {
        prod[i][j] = 0.0;
      }
    }
  }

  if (power == 0) {
    return;  
  }
  // repeated exponentiation
  else {
    long double temp[3][3];
    while (power != 0) {
        if (power % 2 == 1) {
            mat_mul(temp, prod, mat);
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    prod[i][j] = temp[i][j];
                }
            }
        }
        
        power = power >> 1;
        
        if (power == 0) {
            break;
        }

        mat_mul(temp, mat, mat);
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                mat[i][j] = temp[i][j];
            }
        }
    }
  }

  return;
}

// not even any significant speedup using repeated squaring...
int main(int argc, char** argv) {
    long double m1[3][3] = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    long double prod[3][3];

    clock_t start = clock(), diff;
    mat_pow(prod, m1, 500);
    diff = clock() - start;

    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    return 0;
}