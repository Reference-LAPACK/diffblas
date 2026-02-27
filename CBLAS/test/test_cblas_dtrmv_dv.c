/* Test program for cblas_dtrmv forward vector (dv) differentiation */
/* Generated automatically by run_tapenade_cblas.py (same validation as _d and BLAS vector forward) */
/* Mode: dv */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cblas.h"

#ifndef NBDirsMax
#define NBDirsMax 4
#endif
#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE
#define PACKED_SIZE ((MAX_SIZE) * ((MAX_SIZE) + 1) / 2)  /* n*(n+1)/2 for packed storage (match BLAS/test) */

extern void cblas_dtrmv_dv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, CBLAS_INT N, const double *A, double (*Ad)[NBDirsMax], CBLAS_INT lda, double *X, double (*Xd)[NBDirsMax], CBLAS_INT incX, int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    double h = 1.0e-6;  /* Step size for finite differences (match _d test) */
    double atol = 1.0e-5, rtol = 1.0e-5;  /* Pass when abs_error <= atol + rtol*|ad| (same as _d) */
    double max_error = 0.0;  /* max (abs_error/error_bound) over elements (same as _d) */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO Uplo = CblasUpper;
    CBLAS_TRANSPOSE TransA = CblasNoTrans;
    CBLAS_DIAG Diag = CblasNonUnit;
    CBLAS_INT N = TEST_SIZE;
    double A[MAX_SIZE * MAX_SIZE];
    double Ad[MAX_SIZE * MAX_SIZE][NBDirsMax];
    double A_orig[MAX_SIZE * MAX_SIZE];
    double Ad_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT lda = MAX_SIZE;
    double X[MAX_SIZE];
    double Xd[MAX_SIZE][NBDirsMax];
    double X_orig[MAX_SIZE];
    double Xd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incX = 1;
    double X_output[MAX_SIZE];
    double X_ad_output[MAX_SIZE];
    double X_forward[MAX_SIZE], X_backward[MAX_SIZE];

    /* Initialize test data with random numbers (matching _d and Fortran pattern) */
    srand(42);
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) A[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) X[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    /* Initialize derivative seeds (match _d order) */
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Ad[i][idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Xd[i][idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;

    /* Store originals */
    memcpy(A_orig, A, sizeof(A));
    memcpy(Ad_orig, Ad, sizeof(Ad));
    memcpy(X_orig, X, sizeof(X));
    memcpy(Xd_orig, Xd, sizeof(Xd));

    /* Warmup + primal call, save output(s) */
    cblas_dtrmv(
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
    memcpy(X_output, X, sizeof(X));

    /* Restore all inputs and derivative seeds */
    memcpy(A, A_orig, sizeof(A));
    memcpy(Ad, Ad_orig, sizeof(Ad));
    memcpy(X, X_orig, sizeof(X));
    memcpy(Xd, Xd_orig, sizeof(Xd));

    /* Call _dv (implementation uses const void* for alpha/beta in complex, so pass pointers) */
    cblas_dtrmv_dv(
        layout,
        Uplo,
        TransA,
        Diag,
        N,
        A, Ad,
        lda,
        X, Xd,
        incX,
        nbdirs
    );
    memcpy(X_ad_output, X, sizeof(X));

    /* Verify AD primal output matches original (same as _d) */
    {
        double output_diff_max = 0.0;
        for (i = 0; i < MAX_SIZE; i++) {
            double diff = fabs(X_ad_output[i] - X_output[i]);
            if (diff > output_diff_max) output_diff_max = diff;
        }
        if (output_diff_max > 1.0e-10) {
            printf("WARNING: AD function output differs from original (%s): max_diff=%.6e\n", "X", (double)output_diff_max);
        }
    }

    /* Compare results using finite differences (same structure as _d) */
    printf("Testing %s differentiation...\n", "cblas_dtrmv");
    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals (matching _d) */
        memcpy(A, A_orig, sizeof(A));
        memcpy(X, X_orig, sizeof(X));
        /* Forward perturbation: x + h * x_d (same order as _d) */
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) A[j] += h * Ad_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) X[j] += h * Xd_orig[j][idir];
        cblas_dtrmv(
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
        memcpy(X_forward, X, sizeof(X));
        /* Restore primals (matching _d) */
        memcpy(A, A_orig, sizeof(A));
        memcpy(X, X_orig, sizeof(X));
        /* Backward perturbation: x - h * x_d (same order as _d) */
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) A[j] -= h * Ad_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) X[j] -= h * Xd_orig[j][idir];
        cblas_dtrmv(
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
        memcpy(X_backward, X, sizeof(X));
        /* Central diff vs derivative array(s) */
        for (i = 0; i < MAX_SIZE; i++) {
            double fd = (X_forward[i] - X_backward[i]) / (2.0 * h);
            double ad = Xd[i][idir];
            double abs_err = fabs(fd - ad);
            double ad_ref = (fabs(ad) > 1e-10) ? fabs(ad) : 1e-10;
            double bound = atol + rtol * ad_ref;
            if (abs_err > bound) { has_large_errors = 1; }
            double r = abs_err / bound;
            if (r > max_error) max_error = r;
        }
    }
    printf("Maximum error ratio (abs_error/error_bound): %.6e\n", (double)max_error);
    if (has_large_errors) {
        printf("FAIL: Large errors detected in derivatives\n");
        return 1;
    }
    else if (max_error < 0.5) {
        printf("PASS: Derivatives are accurate to machine precision\n");
        return 0;
    }
    else if (max_error < 1.0) {
        printf("PASS: Derivatives are reasonably accurate\n");
        return 0;
    } else {
        printf("WARNING: Derivatives may have significant errors\n");
        return 0;
    }
}

