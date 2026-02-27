/* Test program for cblas_cdotu_sub forward vector (dv) differentiation */
/* Generated automatically by run_tapenade_cblas.py (complex scalar result) */
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

extern void cblas_cdotu_sub_dv(CBLAS_INT N, const void *X, void *Xd, CBLAS_INT incX, const void *Y, void *Yd, CBLAS_INT incY, void *dot, void *dotd, int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    float h = 1.0e-3f;
    float atol = 5.0e-3f, rtol = 5.0e-3f;
    double max_error = 0.0;

    CBLAS_INT N = TEST_SIZE;
    float complex X[MAX_SIZE];
    float complex Xd[MAX_SIZE][NBDirsMax];
    float complex X_orig[MAX_SIZE];
    float complex Xd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incX = 1;
    float complex Y[MAX_SIZE];
    float complex Yd[MAX_SIZE][NBDirsMax];
    float complex Y_orig[MAX_SIZE];
    float complex Yd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incY = 1;
    float complex dot, dot_forward, dot_backward;
    float complex dotd[NBDirsMax];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) X[i] = ((float)rand()/RAND_MAX)*2.0-1.0 + I*((float)rand()/RAND_MAX)*2.0-1.0;
    for (i = 0; i < MAX_SIZE; i++) Y[i] = ((float)rand()/RAND_MAX)*2.0-1.0 + I*((float)rand()/RAND_MAX)*2.0-1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Xd[i][idir] = ((float)rand()/RAND_MAX)*2.0-1.0 + I*((float)rand()/RAND_MAX)*2.0-1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Yd[i][idir] = ((float)rand()/RAND_MAX)*2.0-1.0 + I*((float)rand()/RAND_MAX)*2.0-1.0;

    memcpy(X_orig, X, sizeof(X));
    memcpy(Xd_orig, Xd, sizeof(Xd));
    memcpy(Y_orig, Y, sizeof(Y));
    memcpy(Yd_orig, Yd, sizeof(Yd));

    cblas_cdotu_sub(
        N,
        X,
        incX,
        Y,
        incY,
        &dot
    );

    memcpy(X, X_orig, sizeof(X));
    memcpy(Xd, Xd_orig, sizeof(Xd));
    memcpy(Y, Y_orig, sizeof(Y));
    memcpy(Yd, Yd_orig, sizeof(Yd));

    cblas_cdotu_sub_dv(
        N,
        X, Xd,
        incX,
        Y, Yd,
        incY,
        &dot, dotd,
        nbdirs
    );

    printf("Testing %s differentiation...\n", "cblas_cdotu_sub");
    for (idir = 0; idir < nbdirs; idir++) {
        memcpy(X, X_orig, sizeof(X));
        memcpy(Y, Y_orig, sizeof(Y));
        for (j = 0; j < MAX_SIZE; j++) X[j] += h * Xd_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) Y[j] += h * Yd_orig[j][idir];
        cblas_cdotu_sub(
        N,
        X,
        incX,
        Y,
        incY,
        &dot
    );
        dot_forward = dot;
        memcpy(X, X_orig, sizeof(X));
        memcpy(Y, Y_orig, sizeof(Y));
        for (j = 0; j < MAX_SIZE; j++) X[j] -= h * Xd_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) Y[j] -= h * Yd_orig[j][idir];
        cblas_cdotu_sub(
        N,
        X,
        incX,
        Y,
        incY,
        &dot
    );
        dot_backward = dot;
        float complex fd = (dot_forward - dot_backward) / (2.0 * h);
        float complex ad = dotd[idir];
        double abs_err = cabs(fd - ad);
        double ad_ref = (cabs(ad) > 1e-10) ? cabs(ad) : 1e-10;
        double bound = atol + rtol * ad_ref;
        if (abs_err > bound) { has_large_errors = 1; }
        double r = abs_err / bound;
        if (r > max_error) max_error = r;
    }
    printf("Maximum error ratio (abs_error/error_bound): %.6e\n", max_error);
    if (has_large_errors) { printf("FAIL: Large errors detected in derivatives\n"); return 1; }
    else if (max_error < 0.5) { printf("PASS: Derivatives are accurate to machine precision\n"); return 0; }
    else if (max_error < 2.0) { printf("PASS: Derivatives are reasonably accurate\n"); return 0; }
    else { printf("WARNING: Derivatives may have significant errors\n"); return 0; }
}

