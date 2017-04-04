#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char** argv) {
    int N = 5;
    int P_SIZE = (int) pow(12, N);
    int i, j;

    // initializes all elements to 0

    // long double (*P)[N] = (long double (*)[N]) malloc(P_SIZE*N*sizeof(long double *));

    // for (i = 0; i < 4; i++) {
    //     for (j = 0; j < P_SIZE; j++) {
    //         P[i][j] = 0.0;
    //     }
    // }

    // long double *P[P_SIZE];

    // for (i = 0; i < P_SIZE; i++) {
    //     P[i] = ()
    // }

    long double (*P)[5];

    P = malloc(P_SIZE * sizeof(long double[5]));

    printf("%Ld", P[55][6]);

    return 0;
}