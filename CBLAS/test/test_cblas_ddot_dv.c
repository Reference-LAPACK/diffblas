/* Test program for cblas_ddot forward vector (dv) differentiation */
/* Generated automatically by run_tapenade_cblas.py (scalar result) */
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
#define PACKED_SIZE ((MAX_SIZE) * ((MAX_SIZE) + 1) / 2)

extern void cblas_ddot_dv(CBLAS_INT N, const double *X, double (*Xd)[NBDirsMax], CBLAS_INT incX, const double *Y, double (*Yd)[NBDirsMax], CBLAS_INT incY, double *result, double resultd[NBDirsMax], int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    double h = 1.0e-6;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double max_error = 0.0;

    CBLAS_INT N = TEST_SIZE;
    double X[MAX_SIZE];
    double Xd[MAX_SIZE][NBDirsMax];
    double X_orig[MAX_SIZE];
    double Xd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incX = 1;
    double Y[MAX_SIZE];
    double Yd[MAX_SIZE][NBDirsMax];
    double Y_orig[MAX_SIZE];
    double Yd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incY = 1;
    double result, result_orig;
    double resultd[NBDirsMax];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) X[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) Y[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Xd[i][idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Yd[i][idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;

    memcpy(X_orig, X, sizeof(X));
    memcpy(Xd_orig, Xd, sizeof(Xd));
    memcpy(Y_orig, Y, sizeof(Y));
    memcpy(Yd_orig, Yd, sizeof(Yd));

    result = cblas_ddot(
        N,
        X,
        incX,
        Y,
        incY
    );
    result_orig = result;

    memcpy(X, X_orig, sizeof(X));
    memcpy(Xd, Xd_orig, sizeof(Xd));
    memcpy(Y, Y_orig, sizeof(Y));
    memcpy(Yd, Yd_orig, sizeof(Yd));

    cblas_ddot_dv(
        N,
        X, Xd,
        incX,
        Y, Yd,
        incY,
        &result, resultd,
        nbdirs
    );

    printf("Testing %s differentiation...\n", "cblas_ddot");
    for (idir = 0; idir < nbdirs; idir++) {
        memcpy(X, X_orig, sizeof(X));
        memcpy(Y, Y_orig, sizeof(Y));
        for (j = 0; j < MAX_SIZE; j++) X[j] += h * Xd_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) Y[j] += h * Yd_orig[j][idir];
        double result_forward = cblas_ddot(
        N,
        X,
        incX,
        Y,
        incY
        );
        memcpy(X, X_orig, sizeof(X));
        memcpy(Y, Y_orig, sizeof(Y));
        for (j = 0; j < MAX_SIZE; j++) X[j] -= h * Xd_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) Y[j] -= h * Yd_orig[j][idir];
        double result_backward = cblas_ddot(
        N,
        X,
        incX,
        Y,
        incY
        );
        double fd = (result_forward - result_backward) / (2.0 * h);
        double ad = resultd[idir];
        double abs_err = fabs(fd - ad);
        double ad_ref = (fabs(ad) > 1e-10) ? fabs(ad) : 1e-10;
        double bound = atol + rtol * ad_ref;
        if (abs_err > bound) { has_large_errors = 1; }
        double r = abs_err / bound;
        if (r > max_error) max_error = r;
    }
    printf("Maximum error ratio (abs_error/error_bound): %.6e\n", (double)max_error);
    if (has_large_errors) { printf("FAIL: Large errors detected in derivatives\n"); return 1; }
    else if (max_error < 0.5) { printf("PASS: Derivatives are accurate to machine precision\n"); return 0; }
    else if (max_error < 1.0) { printf("PASS: Derivatives are reasonably accurate\n"); return 0; }
    else { printf("WARNING: Derivatives may have significant errors\n"); return 0; }
}

