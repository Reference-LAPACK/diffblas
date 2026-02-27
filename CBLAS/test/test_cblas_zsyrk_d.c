/* Test program for cblas_zsyrk differentiation */
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
extern void cblas_zsyrk(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE Trans, const CBLAS_INT N, const CBLAS_INT K, const void *alpha, const void *A, const CBLAS_INT lda, const void *beta, void *C, const CBLAS_INT ldc);
/* Differentiated function */
extern void cblas_zsyrk_d(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE Trans, const CBLAS_INT N, const CBLAS_INT K, const void *alpha, void *alpha_d, const void *A, void *A_d, const CBLAS_INT lda, const void *beta, void *beta_d, void *C, void *C_d, const CBLAS_INT ldc);

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
    CBLAS_TRANSPOSE Trans = CblasNoTrans;
    CBLAS_INT N = TEST_SIZE;
    CBLAS_INT K = TEST_SIZE;
    double complex alpha;  /* Will be initialized with random number */
    double complex alpha_orig;  /* Save original value */
    double complex alpha_d;  /* Derivative seed */
    double complex alpha_d_orig;  /* Save derivative seed for finite differences */
    double complex A[MAX_SIZE * MAX_SIZE];
    double complex A_d[MAX_SIZE * MAX_SIZE];  /* Derivative seeds */
    double complex A_d_orig[MAX_SIZE * MAX_SIZE];  /* Save derivative seeds for finite differences */
    double complex A_orig[MAX_SIZE * MAX_SIZE];
    CBLAS_INT lda = MAX_SIZE;
    double complex beta;  /* Will be initialized with random number */
    double complex beta_orig;  /* Save original value */
    double complex beta_d;  /* Derivative seed */
    double complex beta_d_orig;  /* Save derivative seed for finite differences */
    double complex C[MAX_SIZE * MAX_SIZE];
    double complex C_d[MAX_SIZE * MAX_SIZE];  /* Derivative seeds */
    double complex C_d_orig[MAX_SIZE * MAX_SIZE];  /* Save derivative seeds for finite differences */
    double complex C_orig[MAX_SIZE * MAX_SIZE];
    CBLAS_INT ldc = MAX_SIZE;

    /* Initialize test data with random numbers (matching Fortran pattern) */
    srand(42);  /* Seed for reproducibility */
    alpha = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        A[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }
    beta = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        C[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }

    /* Initialize input derivatives to random values (matching Fortran pattern) */
    alpha_d = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        A_d[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }
    beta_d = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        C_d[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + I * (((double)rand() / RAND_MAX) * 2.0 - 1.0);
    }

    /* Store initial derivative values after random initialization (matching Fortran) */
    alpha_d_orig = alpha_d;
    memcpy(A_d_orig, A_d, sizeof(A_d));
    beta_d_orig = beta_d;
    memcpy(C_d_orig, C_d, sizeof(C_d));

    /* Store original values for central difference computation (matching Fortran) */
    alpha_orig = alpha;
    memcpy(A_orig, A, sizeof(A));
    beta_orig = beta;
    memcpy(C_orig, C, sizeof(C));

    /* Call original function */
    cblas_zsyrk(
        layout,
        Uplo,
        Trans,
        N,
        K,
        &alpha,
        A,
        lda,
        &beta,
        C,
        ldc
    );

    /* Save original output */
    double complex C_output[MAX_SIZE * MAX_SIZE];
    memcpy(C_output, C, sizeof(C));

    /* Restore ALL inputs before calling differentiated function */
    /* Note: Derivative seeds were already initialized and saved to _d_orig above */
    alpha = alpha_orig;
    memcpy(A, A_orig, sizeof(A_orig));
    beta = beta_orig;
    memcpy(C, C_orig, sizeof(C_orig));
    /* Restore derivative seeds to ensure they match _d_orig used in finite differences */
    alpha_d = alpha_d_orig;
    memcpy(A_d, A_d_orig, sizeof(A_d_orig));
    beta_d = beta_d_orig;
    memcpy(C_d, C_d_orig, sizeof(C_d_orig));

    /* Call differentiated function with derivative seeds (using _d arrays) */
    cblas_zsyrk_d(
        layout,
        Uplo,
        Trans,
        N,
        K,
        &alpha, &alpha_d,
        A, A_d,
        lda,
        &beta, &beta_d,
        C, C_d,
        ldc
    );

    /* Save AD primal output before FD overwrites C */
    double complex C_ad_output[MAX_SIZE * MAX_SIZE];
    memcpy(C_ad_output, C, sizeof(C));

    /* Compare results using finite differences */
    printf("Testing %s differentiation...\n", "cblas_zsyrk");

    /* Test C derivatives using directional finite differences */
    /* Compute forward and backward perturbations once for all elements */
    double complex C_forward[MAX_SIZE * MAX_SIZE];
    double complex C_backward[MAX_SIZE * MAX_SIZE];

    /* Forward perturbation: x + h * x_d */
    /* Using EXACT same derivative seeds (_d_orig) as in AD call */
    alpha = alpha_orig;
    memcpy(A, A_orig, sizeof(A_orig));
    beta = beta_orig;
    memcpy(C, C_orig, sizeof(C_orig));
    alpha += h * alpha_d_orig;  /* Using EXACT seed from AD call */
    for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) {
        A[j] += h * A_d_orig[j];  /* Using EXACT seed from AD call */
    }
    beta += h * beta_d_orig;  /* Using EXACT seed from AD call */
    for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) {
        C[j] += h * C_d_orig[j];  /* Using EXACT seed from AD call */
    }
    cblas_zsyrk(
        layout,
        Uplo,
        Trans,
        N,
        K,
        &alpha,
        A,
        lda,
        &beta,
        C,
        ldc
    );
    memcpy(C_forward, C, sizeof(C));

    /* Backward perturbation: x - h * x_d */
    /* Using EXACT same derivative seeds (_d_orig) as in AD call */
    alpha = alpha_orig;
    memcpy(A, A_orig, sizeof(A_orig));
    beta = beta_orig;
    memcpy(C, C_orig, sizeof(C_orig));
    alpha -= h * alpha_d_orig;  /* Using EXACT seed from AD call */
    for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) {
        A[j] -= h * A_d_orig[j];  /* Using EXACT seed from AD call */
    }
    beta -= h * beta_d_orig;  /* Using EXACT seed from AD call */
    for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) {
        C[j] -= h * C_d_orig[j];  /* Using EXACT seed from AD call */
    }
    cblas_zsyrk(
        layout,
        Uplo,
        Trans,
        N,
        K,
        &alpha,
        A,
        lda,
        &beta,
        C,
        ldc
    );
    memcpy(C_backward, C, sizeof(C));

    /* Compare AD results with finite differences for each element */
    /* First, verify that AD function produced correct output values (compare saved AD output to original) */
    double output_diff_max = 0.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        double diff = cabs(C_ad_output[i] - C_output[i]);
        if (diff > output_diff_max) output_diff_max = diff;
    }
    if (output_diff_max > 1.0e-10) {
        printf("WARNING: AD function output differs from original: max_diff=%.6e\n", output_diff_max);
    }

    /* Debug: Print first few derivative seeds and AD results */
    printf("Debug: First few derivative seeds and AD results:\n");
    for (i = 0; i < 4; i++) {
        printf("  C_d[%d] = %.6e + %.6e*I, A_d[%d] = %.6e + %.6e*I\n", i, creal(C_d[i]), cimag(C_d[i]), i, creal(A_d_orig[i]), cimag(A_d_orig[i]));
    }
    printf("  alpha_d = %.6e + %.6e*I, beta_d = %.6e + %.6e*I\n", creal(alpha_d_orig), cimag(alpha_d_orig), creal(beta_d_orig), cimag(beta_d_orig));

    /* Check derivatives for output C (all elements) */
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
        double complex fd_derivative = (C_forward[i] - C_backward[i]) / (2.0 * h);
        double complex ad_derivative = C_d[i];

        double ad_mag = cabs(ad_derivative);
        double abs_error = cabs(fd_derivative - ad_derivative);
        double ad_ref = (ad_mag > 1.0e-10) ? ad_mag : 1.0e-10;
        double error_bound = atol + rtol * ad_ref;
        double error_ratio = abs_error / error_bound;  /* > 1 means outside tolerance */
        max_error = (error_ratio > max_error) ? error_ratio : max_error;

        if (error_ratio > 1.0) {
            has_large_errors = 1;
            printf("  Large error in output C[%d]:\n", i);
            printf("    Central diff: %.6e + %.6e*I\n", creal(fd_derivative), cimag(fd_derivative));
            printf("    AD result:   %.6e + %.6e*I\n", creal(ad_derivative), cimag(ad_derivative));
            printf("    Absolute error: %.6e  Error bound: %.6e  Ratio: %.6e\n", abs_error, error_bound, error_ratio);
        }
    }

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