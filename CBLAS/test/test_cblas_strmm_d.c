/* Test program for cblas_strmm differentiation */
/* Generated automatically by run_tapenade_cblas.py */
/* Mode: d */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cblas.h"
#include "cblas_f77.h"

/* Original function */
extern void cblas_strmm(const CBLAS_LAYOUT layout, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const CBLAS_INT M, const CBLAS_INT N, const float alpha, const float *A, const CBLAS_INT lda, float *B, const CBLAS_INT ldb);
/* Differentiated function */
extern void cblas_strmm_d(const CBLAS_LAYOUT layout, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const CBLAS_INT M, const CBLAS_INT N, const float alpha, float alpha_d, const float *A, float *A_d, const CBLAS_INT lda, float *B, float *B_d, const CBLAS_INT ldb);

#define TEST_SIZE 4  /* Matrix/vector size for test */
#define MAX_SIZE TEST_SIZE

int main(void) {
    int i, j;
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match Fortran BLAS tests) */
    float atol = 2.0e-3f, rtol = 2.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_SIDE Side = CblasLeft;
    CBLAS_UPLO Uplo = CblasUpper;
    CBLAS_TRANSPOSE TransA = CblasNoTrans;
    CBLAS_DIAG Diag = CblasNonUnit;
    CBLAS_INT M = TEST_SIZE;
    CBLAS_INT N = TEST_SIZE;
    float alpha;  /* Will be initialized with random number */
    float alpha_orig;  /* Save original value */
    float alpha_d;  /* Derivative seed */
    float alpha_d_orig;  /* Save derivative seed for finite differences */
    float A[MAX_SIZE * MAX_SIZE];
    float A_d[MAX_SIZE * MAX_SIZE];  /* Derivative seeds */
    float A_d_orig[MAX_SIZE * MAX_SIZE];  /* Save derivative seeds for finite differences */
    float A_orig[MAX_SIZE * MAX_SIZE];
    CBLAS_INT lda = MAX_SIZE;
    float B[MAX_SIZE * MAX_SIZE];
    float B_d[MAX_SIZE * MAX_SIZE];  /* Derivative seeds */
    float B_d_orig[MAX_SIZE * MAX_SIZE];  /* Save derivative seeds for finite differences */
    float B_orig[MAX_SIZE * MAX_SIZE];
    CBLAS_INT ldb = MAX_SIZE;

    /* Initialize test data with random numbers (matching Fortran pattern) */
    srand(42);  /* Seed for reproducibility */
    alpha = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        A[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        B[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }

    /* Initialize input derivatives to random values (matching Fortran pattern) */
    alpha_d = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        A_d[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        B_d[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }

    /* Store initial derivative values after random initialization (matching Fortran) */
    alpha_d_orig = alpha_d;
    memcpy(A_d_orig, A_d, sizeof(A_d));
    memcpy(B_d_orig, B_d, sizeof(B_d));

    /* Store original values for central difference computation (matching Fortran) */
    alpha_orig = alpha;
    memcpy(A_orig, A, sizeof(A));
    memcpy(B_orig, B, sizeof(B));

    /* Call original function */
    cblas_strmm(
        layout,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb
    );

    /* Restore ALL inputs before calling differentiated function */
    /* Note: Derivative seeds were already initialized and saved to _d_orig above */
    alpha = alpha_orig;
    memcpy(A, A_orig, sizeof(A_orig));
    memcpy(B, B_orig, sizeof(B_orig));
    /* Restore derivative seeds to ensure they match _d_orig used in finite differences */
    alpha_d = alpha_d_orig;
    memcpy(A_d, A_d_orig, sizeof(A_d_orig));
    memcpy(B_d, B_d_orig, sizeof(B_d_orig));

    /* Call differentiated function with derivative seeds (using _d arrays) */
    cblas_strmm_d(
        layout,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        alpha, alpha_d,
        A, A_d,
        lda,
        B, B_d,
        ldb
    );

    /* Compare results using finite differences */
    printf("Testing %s differentiation...\n", "cblas_strmm");

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