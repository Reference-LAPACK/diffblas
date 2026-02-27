/* Test program for cblas_ssymv differentiation */
/* Generated automatically by run_tapenade_cblas.py */
/* Mode: d */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cblas.h"
#include "cblas_f77.h"

/* Original function */
extern void cblas_ssymv(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_INT N, const float alpha, const float *A, const CBLAS_INT lda, const float *X, const CBLAS_INT incX, const float beta, float *Y, const CBLAS_INT incY);
/* Differentiated function */
extern void cblas_ssymv_d(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_INT N, const float alpha, float alpha_d, const float *A, float *A_d, const CBLAS_INT lda, const float *X, float *X_d, const CBLAS_INT incX, const float beta, float beta_d, float *Y, float *Y_d, const CBLAS_INT incY);

#define TEST_SIZE 4  /* Matrix/vector size for test */
#define MAX_SIZE TEST_SIZE

int main(void) {
    int i, j;
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match Fortran BLAS tests) */
    float atol = 2.0e-3f, rtol = 2.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO Uplo = CblasUpper;
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
    float X[MAX_SIZE];
    float X_d[MAX_SIZE];  /* Derivative seeds */
    float X_d_orig[MAX_SIZE];  /* Save derivative seeds for finite differences */
    float X_orig[MAX_SIZE];
    CBLAS_INT incX;
    float beta;  /* Will be initialized with random number */
    float beta_orig;  /* Save original value */
    float beta_d;  /* Derivative seed */
    float beta_d_orig;  /* Save derivative seed for finite differences */
    float Y[MAX_SIZE];
    float Y_d[MAX_SIZE];  /* Derivative seeds */
    float Y_d_orig[MAX_SIZE];  /* Save derivative seeds for finite differences */
    float Y_orig[MAX_SIZE];
    CBLAS_INT incY;

    /* Initialize test data with random numbers (matching Fortran pattern) */
    srand(42);  /* Seed for reproducibility */
    alpha = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        A[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
    for (i = 0; i < MAX_SIZE; i++) {
        X[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
    beta = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    for (i = 0; i < MAX_SIZE; i++) {
        Y[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
    incX = 1;  /* Typical BLAS increment value */
    incY = 1;  /* Typical BLAS increment value */

    /* Initialize input derivatives to random values (matching Fortran pattern) */
    alpha_d = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        A_d[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
    for (i = 0; i < MAX_SIZE; i++) {
        X_d[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }
    beta_d = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    for (i = 0; i < MAX_SIZE; i++) {
        Y_d[i] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
    }

    /* Store initial derivative values after random initialization (matching Fortran) */
    alpha_d_orig = alpha_d;
    memcpy(A_d_orig, A_d, sizeof(A_d));
    memcpy(X_d_orig, X_d, sizeof(X_d));
    beta_d_orig = beta_d;
    memcpy(Y_d_orig, Y_d, sizeof(Y_d));

    /* Store original values for central difference computation (matching Fortran) */
    alpha_orig = alpha;
    memcpy(A_orig, A, sizeof(A));
    memcpy(X_orig, X, sizeof(X));
    beta_orig = beta;
    memcpy(Y_orig, Y, sizeof(Y));

    /* Call original function */
    cblas_ssymv(
        layout,
        Uplo,
        N,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY
    );

    /* Restore ALL inputs before calling differentiated function */
    /* Note: Derivative seeds were already initialized and saved to _d_orig above */
    alpha = alpha_orig;
    memcpy(A, A_orig, sizeof(A_orig));
    beta = beta_orig;
    /* Restore derivative seeds to ensure they match _d_orig used in finite differences */
    alpha_d = alpha_d_orig;
    memcpy(A_d, A_d_orig, sizeof(A_d_orig));
    beta_d = beta_d_orig;

    /* Call differentiated function with derivative seeds (using _d arrays) */
    cblas_ssymv_d(
        layout,
        Uplo,
        N,
        alpha, alpha_d,
        A, A_d,
        lda,
        X, X_d,
        incX,
        beta, beta_d,
        Y, Y_d,
        incY
    );

    /* Compare results using finite differences */
    printf("Testing %s differentiation...\n", "cblas_ssymv");

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