/* Test program for cblas_ztrmv differentiation */
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
extern void cblas_ztrmv(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const CBLAS_INT N, const void *A, const CBLAS_INT lda, void *X, const CBLAS_INT incX);
/* Differentiated function */
extern void cblas_ztrmv_d(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const CBLAS_INT N, const void *A, void *A_d, const CBLAS_INT lda, void *X, void *X_d, const CBLAS_INT incX);

#define TEST_SIZE 4  /* Matrix/vector size for test */
#define MAX_SIZE TEST_SIZE

int main(void) {
    int i, j;
    int has_large_errors = 0;
    double h = 1.0e-6;  /* Step size for finite differences (match Fortran BLAS tests) */
    double atol = 1.0e-5, rtol = 1.0e-5;  /* Pass when abs_error <= atol + rtol*|ad| */
    double max_error = 0.0;  /* max (abs_error/error_bound) over elements */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO Uplo = CblasUpper;
    CBLAS_TRANSPOSE TransA = CblasNoTrans;
    CBLAS_DIAG Diag = CblasNonUnit;
    CBLAS_INT N = TEST_SIZE;
    double complex A[MAX_SIZE * MAX_SIZE];
    double complex A_d[MAX_SIZE * MAX_SIZE];  /* Derivative seeds */
    double complex A_d_orig[MAX_SIZE * MAX_SIZE];  /* Save derivative seeds for finite differences */
    double complex A_orig[MAX_SIZE * MAX_SIZE];
    CBLAS_INT lda = MAX_SIZE;
    double complex X[MAX_SIZE];
    double complex X_d[MAX_SIZE];  /* Derivative seeds */
    double complex X_d_orig[MAX_SIZE];  /* Save derivative seeds for finite differences */
    double complex X_orig[MAX_SIZE];
    CBLAS_INT incX;

    /* Initialize test data with random numbers (matching Fortran pattern) */
    srand(42);  /* Seed for reproducibility */
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        A[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }
    for (i = 0; i < MAX_SIZE; i++) {
        X[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }
    incX = 1;  /* Typical BLAS increment value */

    /* Initialize input derivatives to random values (matching Fortran pattern) */
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        A_d[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }
    for (i = 0; i < MAX_SIZE; i++) {
        X_d[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }

    /* Store initial derivative values after random initialization (matching Fortran) */
    memcpy(A_d_orig, A_d, sizeof(A_d));
    memcpy(X_d_orig, X_d, sizeof(X_d));

    /* Store original values for central difference computation (matching Fortran) */
    memcpy(A_orig, A, sizeof(A));
    memcpy(X_orig, X, sizeof(X));

    /* Call original function */
    cblas_ztrmv(
        layout,
        Uplo,
        TransA,
        Diag,
        N,
        A,
        lda,
        X,
        incX
    );

    /* Restore ALL inputs before calling differentiated function */
    /* Note: Derivative seeds were already initialized and saved to _d_orig above */
    memcpy(A, A_orig, sizeof(A_orig));
    /* Restore derivative seeds to ensure they match _d_orig used in finite differences */
    memcpy(A_d, A_d_orig, sizeof(A_d_orig));

    /* Call differentiated function with derivative seeds (using _d arrays) */
    cblas_ztrmv_d(
        layout,
        Uplo,
        TransA,
        Diag,
        N,
        A, A_d,
        lda,
        X, X_d,
        incX
    );

    /* Compare results using finite differences */
    printf("Testing %s differentiation...\n", "cblas_ztrmv");

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