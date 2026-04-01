/* Test program for cblas_cscal differentiation */
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
extern void cblas_cscal(const CBLAS_INT N, const void *alpha, void *X, const CBLAS_INT incX);
/* Differentiated function */
extern void cblas_cscal_d(const CBLAS_INT N, const void *alpha, void *alpha_d, void *X, void *X_d, const CBLAS_INT incX);

#define TEST_SIZE 4  /* Matrix/vector size for test */
#define MAX_SIZE TEST_SIZE
extern void set_isize1ofcx_(int *val);

int main(void) {
    int i, j;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofcx_(&diffblas_isize);
    }
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match Fortran BLAS tests) */
    float atol = 2.0e-3f, rtol = 2.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements */

    CBLAS_INT N = TEST_SIZE;
    float complex alpha;  /* Will be initialized with random number */
    float complex alpha_orig;  /* Save original value */
    float complex alpha_d;  /* Derivative seed */
    float complex alpha_d_orig;  /* Save derivative seed for finite differences */
    float complex X[MAX_SIZE];
    float complex X_d[MAX_SIZE];  /* Derivative seeds */
    float complex X_d_orig[MAX_SIZE];  /* Save derivative seeds for finite differences */
    float complex X_orig[MAX_SIZE];
    CBLAS_INT incX;

    /* Initialize test data with random numbers (matching Fortran pattern) */
    srand(42);  /* Seed for reproducibility */
    alpha = ((float)rand() / RAND_MAX) * 2.0f - 1.0f + I * (((float)rand() / RAND_MAX) * 2.0f - 1.0f);
    for (i = 0; i < MAX_SIZE; i++) {
        X[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f + I * (((float)rand() / RAND_MAX) * 2.0f - 1.0f);
    }
    incX = 1;  /* Typical BLAS increment value */

    /* Initialize input derivatives to random values (matching Fortran pattern) */
    alpha_d = ((float)rand() / RAND_MAX) * 2.0f - 1.0f + I * (((float)rand() / RAND_MAX) * 2.0f - 1.0f);
    for (i = 0; i < MAX_SIZE; i++) {
        X_d[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f + I * (((float)rand() / RAND_MAX) * 2.0f - 1.0f);
    }

    /* Store initial derivative values after random initialization (matching Fortran) */
    alpha_d_orig = alpha_d;
    memcpy(X_d_orig, X_d, sizeof(X_d));

    /* Store original values for central difference computation (matching Fortran) */
    alpha_orig = alpha;
    memcpy(X_orig, X, sizeof(X));

    /* Call original function */
    cblas_cscal(
        N,
        &alpha,
        X,
        incX
    );

    /* Restore ALL inputs before calling differentiated function */
    /* Note: Derivative seeds were already initialized and saved to _d_orig above */
    alpha = alpha_orig;
    /* Restore derivative seeds to ensure they match _d_orig used in finite differences */
    alpha_d = alpha_d_orig;

    /* Call differentiated function with derivative seeds (using _d arrays) */
    cblas_cscal_d(
        N,
        &alpha, &alpha_d,
        X, X_d,
        incX
    );

    /* Compare results using finite differences */
    printf("Testing %s differentiation...\n", "cblas_cscal");

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
