/* Test program for cblas_zswap forward vector (dv) differentiation */
/* Generated automatically by run_tapenade_cblas.py (same validation as _d and BLAS vector forward) */
/* Mode: dv */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "cblas.h"

#ifndef NBDirsMax
#define NBDirsMax 4
#endif
#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE
#define PACKED_SIZE ((MAX_SIZE) * ((MAX_SIZE) + 1) / 2)  /* n*(n+1)/2 for packed storage (match BLAS/test) */

extern void cblas_zswap_dv(CBLAS_INT N, void *X, void *Xd, CBLAS_INT incX, void *Y, void *Yd, CBLAS_INT incY, int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    double h = 1.0e-6;  /* Step size for finite differences (match _d test) */
    double atol = 1.0e-5, rtol = 1.0e-5;  /* Pass when abs_error <= atol + rtol*|ad| (same as _d) */
    double max_error = 0.0;  /* max (abs_error/error_bound) over elements (same as _d) */

    CBLAS_INT N = TEST_SIZE;
    double complex X[MAX_SIZE];
    double complex Xd[MAX_SIZE][NBDirsMax];
    double complex X_orig[MAX_SIZE];
    double complex Xd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incX = 1;
    double complex Y[MAX_SIZE];
    double complex Yd[MAX_SIZE][NBDirsMax];
    double complex Y_orig[MAX_SIZE];
    double complex Yd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incY = 1;
    double complex X_output[MAX_SIZE];
    double complex X_ad_output[MAX_SIZE];
    double complex X_forward[MAX_SIZE], X_backward[MAX_SIZE];
    double complex Y_output[MAX_SIZE];
    double complex Y_ad_output[MAX_SIZE];
    double complex Y_forward[MAX_SIZE], Y_backward[MAX_SIZE];

    /* Initialize test data with random numbers (matching _d and Fortran pattern) */
    srand(42);
    for (i = 0; i < MAX_SIZE; i++) X[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) Y[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    /* Initialize derivative seeds (match _d order) */
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Xd[i][idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Yd[i][idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;

    /* Store originals */
    memcpy(X_orig, X, sizeof(X));
    memcpy(Xd_orig, Xd, sizeof(Xd));
    memcpy(Y_orig, Y, sizeof(Y));
    memcpy(Yd_orig, Yd, sizeof(Yd));

    /* Warmup + primal call, save output(s) */
    cblas_zswap(
        N,
        X,
        incX,
        Y,
        incY
    );
    memcpy(X_output, X, sizeof(X));
    memcpy(Y_output, Y, sizeof(Y));

    /* Restore all inputs and derivative seeds */
    memcpy(X, X_orig, sizeof(X));
    memcpy(Xd, Xd_orig, sizeof(Xd));
    memcpy(Y, Y_orig, sizeof(Y));
    memcpy(Yd, Yd_orig, sizeof(Yd));

    /* Call _dv (implementation uses const void* for alpha/beta in complex, so pass pointers) */
    cblas_zswap_dv(
        N,
        X, Xd,
        incX,
        Y, Yd,
        incY,
        nbdirs
    );
    memcpy(X_ad_output, X, sizeof(X));
    memcpy(Y_ad_output, Y, sizeof(Y));

    /* Verify AD primal output matches original (same as _d) */
    {
        double output_diff_max = 0.0;
        for (i = 0; i < MAX_SIZE; i++) {
            double diff = cabs(X_ad_output[i] - X_output[i]);
            if (diff > output_diff_max) output_diff_max = diff;
        }
        if (output_diff_max > 1.0e-10) {
            printf("WARNING: AD function output differs from original (%s): max_diff=%.6e\n", "X", (double)output_diff_max);
        }
    }
    {
        double output_diff_max = 0.0;
        for (i = 0; i < MAX_SIZE; i++) {
            double diff = cabs(Y_ad_output[i] - Y_output[i]);
            if (diff > output_diff_max) output_diff_max = diff;
        }
        if (output_diff_max > 1.0e-10) {
            printf("WARNING: AD function output differs from original (%s): max_diff=%.6e\n", "Y", (double)output_diff_max);
        }
    }

    /* Compare results using finite differences (same structure as _d) */
    printf("Testing %s differentiation...\n", "cblas_zswap");
    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals (matching _d) */
        memcpy(X, X_orig, sizeof(X));
        memcpy(Y, Y_orig, sizeof(Y));
        /* Forward perturbation: x + h * x_d (same order as _d) */
        for (j = 0; j < MAX_SIZE; j++) X[j] += h * Xd_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) Y[j] += h * Yd_orig[j][idir];
        cblas_zswap(
        N,
        X,
        incX,
        Y,
        incY
        );
        memcpy(X_forward, X, sizeof(X));
        memcpy(Y_forward, Y, sizeof(Y));
        /* Restore primals (matching _d) */
        memcpy(X, X_orig, sizeof(X));
        memcpy(Y, Y_orig, sizeof(Y));
        /* Backward perturbation: x - h * x_d (same order as _d) */
        for (j = 0; j < MAX_SIZE; j++) X[j] -= h * Xd_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) Y[j] -= h * Yd_orig[j][idir];
        cblas_zswap(
        N,
        X,
        incX,
        Y,
        incY
        );
        memcpy(X_backward, X, sizeof(X));
        memcpy(Y_backward, Y, sizeof(Y));
        /* Central diff vs derivative array(s) */
        for (i = 0; i < MAX_SIZE; i++) {
            double fd = (X_forward[i] - X_backward[i]) / (2.0 * h);
            double ad = Xd[i][idir];
            double abs_err = cabs(fd - ad);
            double ad_ref = (cabs(ad) > 1e-10) ? cabs(ad) : 1e-10;
            double bound = atol + rtol * ad_ref;
            if (abs_err > bound) { has_large_errors = 1; }
            double r = abs_err / bound;
            if (r > max_error) max_error = r;
        }
        for (i = 0; i < MAX_SIZE; i++) {
            double fd = (Y_forward[i] - Y_backward[i]) / (2.0 * h);
            double ad = Yd[i][idir];
            double abs_err = cabs(fd - ad);
            double ad_ref = (cabs(ad) > 1e-10) ? cabs(ad) : 1e-10;
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

