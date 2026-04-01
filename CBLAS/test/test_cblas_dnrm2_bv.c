/* Test program for cblas_dnrm2 vector reverse (bv) differentiation (scalar result) */
/* Generated automatically by run_tapenade_cblas.py - same validation as _dv scalar */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cblas.h"
#include "cblas_f77.h"
#include "cblas_bv.h"

#ifndef NBDirsMax
#define NBDirsMax 4
#endif
#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE

extern double cblas_dnrm2(CBLAS_INT N, const double *X, CBLAS_INT incX);
extern void cblas_dnrm2_bv(CBLAS_INT N, const double *X, double (*X_b)[NBDirsMax], CBLAS_INT incX, double result_b[NBDirsMax], int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    double h = 1.0e-6;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double max_error = 0.0;
    double vjp_fd, vjp_ad;

    CBLAS_INT N = TEST_SIZE;
    double X[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    double X_b[MAX_SIZE][NBDirsMax];
    CBLAS_INT incX = 1;
    double result_b[NBDirsMax];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) X[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) X_dir[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;

    memcpy(X_orig, X, sizeof(X));

    for (i = 0; i < MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) X_b[i][j] = 0.0;
    for (j = 0; j < NBDirsMax; j++) result_b[j] = 1.0;  /* seed cotangent for scalar result */

    cblas_dnrm2_bv(N, X, X_b, incX, result_b, nbdirs);

    for (idir = 0; idir < nbdirs; idir++) {
        memcpy(X, X_orig, sizeof(X));
        for (i = 0; i < MAX_SIZE; i++) X_dir[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        /* Forward */
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
        double result_forward = cblas_dnrm2(
        N,
        X,
        incX
        );
        memcpy(X, X_orig, sizeof(X));
        /* Backward */
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
        double result_backward = cblas_dnrm2(
        N,
        X,
        incX
        );
        vjp_fd = (result_forward - result_backward) / (2.0 * h);
        vjp_ad = 0.0;
        for (i = 0; i < MAX_SIZE; i++) vjp_ad += X_dir[i] * X_b[i][idir];
        {
            double abs_err = fabs(vjp_fd - vjp_ad);
            double ref = (fabs(vjp_ad) > 1e-10) ? fabs(vjp_ad) : 1e-10;
            double bound = atol + rtol * ref;
            if (abs_err > bound) has_large_errors = 1;
            { double r = abs_err / bound; if (r > max_error) max_error = r; }
        }
    }
    printf("Maximum error ratio: %.6e\n", (double)max_error);
    if (has_large_errors) { printf("FAIL: Large errors in derivatives\n"); return 1; }
    if (max_error < 0.5) { printf("PASS: Derivatives accurate to machine precision\n"); return 0; }
    if (max_error < 1.0) { printf("PASS: Derivatives reasonably accurate\n"); return 0; }
    printf("WARNING: Derivatives may have significant errors\n"); return 0;
}

