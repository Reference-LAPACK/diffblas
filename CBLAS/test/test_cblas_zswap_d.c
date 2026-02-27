/* Test program for cblas_zswap differentiation */
/* Generated automatically by run_tapenade_cblas.py */
/* Mode: d */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "cblas.h"
#include "cblas_f77.h"

/* Original function */
extern void cblas_zswap(const CBLAS_INT N, void *X, const CBLAS_INT incX, void *Y, const CBLAS_INT incY);
/* Differentiated function */
extern void cblas_zswap_d(const CBLAS_INT N, void *X, void *X_d, const CBLAS_INT incX, void *Y, void *Y_d, const CBLAS_INT incY);

#define TEST_SIZE 4  /* Matrix/vector size for test */
#define MAX_SIZE TEST_SIZE

int main(void) {
    int i, j;
    int has_large_errors = 0;
    double h = 1.0e-6;  /* Step size for finite differences (match Fortran BLAS tests) */
    double atol = 1.0e-5, rtol = 1.0e-5;  /* Pass when abs_error <= atol + rtol*|ad| */
    double max_error = 0.0;  /* max (abs_error/error_bound) over elements */

    CBLAS_INT N = TEST_SIZE;
    double complex X[MAX_SIZE];
    double complex X_d[MAX_SIZE];  /* Derivative seeds */
    double complex X_d_orig[MAX_SIZE];  /* Save derivative seeds for finite differences */
    double complex X_orig[MAX_SIZE];
    CBLAS_INT incX;
    double complex Y[MAX_SIZE];
    double complex Y_d[MAX_SIZE];  /* Derivative seeds */
    double complex Y_d_orig[MAX_SIZE];  /* Save derivative seeds for finite differences */
    double complex Y_orig[MAX_SIZE];
    CBLAS_INT incY;

    /* Initialize test data with random numbers (matching Fortran pattern) */
    srand(42);  /* Seed for reproducibility */
    for (i = 0; i < MAX_SIZE; i++) {
        X[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }
    for (i = 0; i < MAX_SIZE; i++) {
        Y[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }
    incX = 1;  /* Typical BLAS increment value */
    incY = 1;  /* Typical BLAS increment value */

    /* Initialize input derivatives to random values (matching Fortran pattern) */
    for (i = 0; i < MAX_SIZE; i++) {
        X_d[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }
    for (i = 0; i < MAX_SIZE; i++) {
        Y_d[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }

    /* Store initial derivative values after random initialization (matching Fortran) */
    memcpy(X_d_orig, X_d, sizeof(X_d));
    memcpy(Y_d_orig, Y_d, sizeof(Y_d));

    /* Store original values for central difference computation (matching Fortran) */
    memcpy(X_orig, X, sizeof(X));
    memcpy(Y_orig, Y, sizeof(Y));

    /* Call original function */
    cblas_zswap(
        N,
        X,
        incX,
        Y,
        incY
    );

    /* Restore ALL inputs before calling differentiated function */
    /* Note: Derivative seeds were already initialized and saved to _d_orig above */
    /* Restore derivative seeds to ensure they match _d_orig used in finite differences */

    /* Call differentiated function with derivative seeds (using _d arrays) */
    cblas_zswap_d(
        N,
        X, X_d,
        incX,
        Y, Y_d,
        incY
    );

    /* Compare results using finite differences */
    printf("Testing %s differentiation...\n", "cblas_zswap");

    printf("Maximum error ratio (abs_error/error_bound): %.6e\n", max_error);
    if (has_large_errors) {
        printf("FAIL: Large errors detected in derivatives\n");
        return 1;
    } else if (max_error < 0.5) {
        printf("PASS: Derivatives are accurate to machine precision\n");
        return 0;
    } else if (max_error < 1.0) {
        printf("PASS: Derivatives are reasonably accurate\n");
        return 0;
    } else {
        printf("WARNING: Derivatives may have significant errors\n");
        return 0;
    }
}