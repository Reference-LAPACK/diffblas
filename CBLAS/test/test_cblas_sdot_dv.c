/* Test program for cblas_sdot forward vector (dv) differentiation */
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

extern void cblas_sdot_dv(CBLAS_INT N, const float *X, float (*Xd)[NBDirsMax], CBLAS_INT incX, const float *Y, float (*Yd)[NBDirsMax], CBLAS_INT incY, float *result, float resultd[NBDirsMax], int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    float h = 1.0e-3f;
    float atol = 5.0e-3f, rtol = 5.0e-3f;
    float max_error = 0.0f;

    CBLAS_INT N = TEST_SIZE;
    float X[MAX_SIZE];
    float Xd[MAX_SIZE][NBDirsMax];
    float X_orig[MAX_SIZE];
    float Xd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incX = 1;
    float Y[MAX_SIZE];
    float Yd[MAX_SIZE][NBDirsMax];
    float Y_orig[MAX_SIZE];
    float Yd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incY = 1;
    float result, result_orig;
    float resultd[NBDirsMax];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) X[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) Y[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Xd[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Yd[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;

    memcpy(X_orig, X, sizeof(X));
    memcpy(Xd_orig, Xd, sizeof(Xd));
    memcpy(Y_orig, Y, sizeof(Y));
    memcpy(Yd_orig, Yd, sizeof(Yd));

    result = cblas_sdot(
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

    cblas_sdot_dv(
        N,
        X, Xd,
        incX,
        Y, Yd,
        incY,
        &result, resultd,
        nbdirs
    );

    printf("Testing %s differentiation...\n", "cblas_sdot");
    for (idir = 0; idir < nbdirs; idir++) {
        memcpy(X, X_orig, sizeof(X));
        memcpy(Y, Y_orig, sizeof(Y));
        for (j = 0; j < MAX_SIZE; j++) X[j] += h * Xd_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) Y[j] += h * Yd_orig[j][idir];
        float result_forward = cblas_sdot(
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
        float result_backward = cblas_sdot(
        N,
        X,
        incX,
        Y,
        incY
        );
        float fd = (result_forward - result_backward) / (2.0 * h);
        float ad = resultd[idir];
        float abs_err = fabs(fd - ad);
        float ad_ref = (fabs(ad) > 1e-10) ? fabs(ad) : 1e-10;
        float bound = atol + rtol * ad_ref;
        if (abs_err > bound) { has_large_errors = 1; }
        float r = abs_err / bound;
        if (r > max_error) max_error = r;
    }
    printf("Maximum error ratio (abs_error/error_bound): %.6e\n", (double)max_error);
    if (has_large_errors) { printf("FAIL: Large errors detected in derivatives\n"); return 1; }
    else if (max_error < 0.5f) { printf("PASS: Derivatives are accurate to machine precision\n"); return 0; }
    else if (max_error < 2.0f) { printf("PASS: Derivatives are reasonably accurate\n"); return 0; }
    else { printf("WARNING: Derivatives may have significant errors\n"); return 0; }
}

