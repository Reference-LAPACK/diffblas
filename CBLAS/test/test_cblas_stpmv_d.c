/* Test program for cblas_stpmv differentiation */
/* Generated automatically by run_tapenade_cblas.py */
/* Mode: d */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cblas.h"
#include "cblas_f77.h"

/* Original function */
extern void cblas_stpmv(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const CBLAS_INT N, const float *Ap, float *X, const CBLAS_INT incX);
/* Differentiated function */
extern void cblas_stpmv_d(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const CBLAS_INT N, const float *Ap, float *Ap_d, float *X, float *X_d, const CBLAS_INT incX);

#define TEST_SIZE 4  /* Matrix/vector size for test */
#define MAX_SIZE TEST_SIZE
#define PACKED_SIZE ((MAX_SIZE) * ((MAX_SIZE) + 1) / 2)  /* packed symmetric/triangular */

int main(void) {
    int i, j;
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match Fortran BLAS tests) */
    float atol = 2.0e-3f, rtol = 2.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO Uplo = CblasUpper;
    CBLAS_TRANSPOSE TransA = CblasNoTrans;
    CBLAS_DIAG Diag = CblasNonUnit;
    CBLAS_INT N = TEST_SIZE;
    float Ap[PACKED_SIZE];
    float Ap_d[PACKED_SIZE];  /* Derivative seeds */
    float Ap_d_orig[PACKED_SIZE];
    float Ap_orig[PACKED_SIZE];
    float X[MAX_SIZE];
    float X_d[MAX_SIZE];  /* Derivative seeds */
    float X_d_orig[MAX_SIZE];  /* Save derivative seeds for finite differences */
    float X_orig[MAX_SIZE];
    CBLAS_INT incX;

    /* Initialize test data with random numbers (matching Fortran pattern) */
    srand(42);  /* Seed for reproducibility */
    for (j = 0; j < MAX_SIZE; j++) {
        for (i = 0; i <= j; i++) {
            Ap[j * (j + 1) / 2 + i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
        }
    }
    for (i = 0; i < MAX_SIZE; i++) {
        X[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
    incX = 1;  /* Typical BLAS increment value */

    /* Initialize input derivatives to random values (matching Fortran pattern) */
    for (i = 0; i < PACKED_SIZE; i++) {
        Ap_d[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
    for (i = 0; i < MAX_SIZE; i++) {
        X_d[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }

    /* Store initial derivative values after random initialization (matching Fortran) */
    memcpy(Ap_d_orig, Ap_d, sizeof(Ap_d));
    memcpy(X_d_orig, X_d, sizeof(X_d));

    /* Store original values for central difference computation (matching Fortran) */
    memcpy(Ap_orig, Ap, sizeof(Ap));
    memcpy(X_orig, X, sizeof(X));

    /* Call original function */
    cblas_stpmv(
        layout,
        Uplo,
        TransA,
        Diag,
        N,
        Ap,
        X,
        incX
    );

    /* Restore ALL inputs before calling differentiated function */
    /* Note: Derivative seeds were already initialized and saved to _d_orig above */
    memcpy(Ap, Ap_orig, sizeof(Ap_orig));
    /* Restore derivative seeds to ensure they match _d_orig used in finite differences */
    memcpy(Ap_d, Ap_d_orig, sizeof(Ap_d_orig));

    /* Call differentiated function with derivative seeds (using _d arrays) */
    cblas_stpmv_d(
        layout,
        Uplo,
        TransA,
        Diag,
        N,
        Ap, Ap_d,
        X, X_d,
        incX
    );

    /* Compare results using finite differences */
    printf("Testing %s differentiation...\n", "cblas_stpmv");

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