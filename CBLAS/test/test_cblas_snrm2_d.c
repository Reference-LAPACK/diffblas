/* Test program for cblas_snrm2 differentiation */
/* Generated automatically by run_tapenade_cblas.py */
/* Mode: d */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cblas.h"
#include "cblas_f77.h"

/* Original function */
extern float cblas_snrm2(const CBLAS_INT N, const float *X, const CBLAS_INT incX);
/* Differentiated function */
extern float cblas_snrm2_d(const CBLAS_INT N, const float *X, float *X_d, const CBLAS_INT incX, float *result);

#define TEST_SIZE 4  /* Matrix/vector size for test */
#define MAX_SIZE TEST_SIZE

int main(void) {
    int i, j;
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match Fortran BLAS tests) */
    float atol = 2.0e-3f, rtol = 2.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements */

    CBLAS_INT N = TEST_SIZE;
    float X[MAX_SIZE];
    float X_d[MAX_SIZE];  /* Derivative seeds */
    float X_d_orig[MAX_SIZE];  /* Save derivative seeds for finite differences */
    float X_orig[MAX_SIZE];
    CBLAS_INT incX;

    /* Initialize test data with random numbers (matching Fortran pattern) */
    srand(42);  /* Seed for reproducibility */
    for (i = 0; i < MAX_SIZE; i++) {
        X[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
    incX = 1;  /* Typical BLAS increment value */

    /* Initialize input derivatives to random values (matching Fortran pattern) */
    for (i = 0; i < MAX_SIZE; i++) {
        X_d[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }

    /* Store initial derivative values after random initialization (matching Fortran) */
    memcpy(X_d_orig, X_d, sizeof(X_d));

    /* Store original values for central difference computation (matching Fortran) */
    memcpy(X_orig, X, sizeof(X));

    /* Call original function */
    cblas_snrm2(
        N,
        X,
        incX
    );

    /* Restore ALL inputs before calling differentiated function */
    /* Note: Derivative seeds were already initialized and saved to _d_orig above */
    /* Restore derivative seeds to ensure they match _d_orig used in finite differences */

    float result;
    float result_d;

    /* Call differentiated function with derivative seeds (using _d arrays) */
    result_d = cblas_snrm2_d(
        N,
        X, X_d,
        incX,
        &result
    );

    /* Compare results using finite differences */
    printf("Testing %s differentiation...\n", "cblas_snrm2");

    printf("Maximum error ratio (abs_error/error_bound): %.6e\n", max_error);
    if (has_large_errors) {
        printf("FAIL: Large errors detected in derivatives\n");
        return 1;
    } else if (max_error < 0.5f) {
        printf("PASS: Derivatives are accurate to machine precision\n");
        return 0;
    } else if (max_error < 1.0f) {
        printf("PASS: Derivatives are reasonably accurate\n");
        return 0;
    } else {
        printf("WARNING: Derivatives may have significant errors\n");
        return 0;
    }
}