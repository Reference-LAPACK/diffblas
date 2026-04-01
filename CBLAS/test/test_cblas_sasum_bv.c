/* Test program for cblas_sasum vector reverse (bv) differentiation (scalar result) */
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

extern float cblas_sasum(CBLAS_INT N, const float *X, CBLAS_INT incX);
extern void cblas_sasum_bv(CBLAS_INT N, const float *X, float (*X_b)[NBDirsMax], CBLAS_INT incX, float result_b[NBDirsMax], int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    float h = 1.0e-3f;
    float atol = 2.0e-3f, rtol = 2.0e-3f;
    float max_error = 0.0f;
    float vjp_fd, vjp_ad;

    CBLAS_INT N = TEST_SIZE;
    float X[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    float X_b[MAX_SIZE][NBDirsMax];
    CBLAS_INT incX = 1;
    float result_b[NBDirsMax];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) X[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) X_dir[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;

    memcpy(X_orig, X, sizeof(X));

    for (i = 0; i < MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) X_b[i][j] = 0.0f;
    for (j = 0; j < NBDirsMax; j++) result_b[j] = 1.0f;  /* seed cotangent for scalar result */

    cblas_sasum_bv(N, X, X_b, incX, result_b, nbdirs);

    for (idir = 0; idir < nbdirs; idir++) {
        memcpy(X, X_orig, sizeof(X));
        for (i = 0; i < MAX_SIZE; i++) X_dir[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
        /* Forward */
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
        float result_forward = cblas_sasum(
        N,
        X,
        incX
        );
        memcpy(X, X_orig, sizeof(X));
        /* Backward */
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
        float result_backward = cblas_sasum(
        N,
        X,
        incX
        );
        vjp_fd = (result_forward - result_backward) / (2.0 * h);
        vjp_ad = 0.0f;
        for (i = 0; i < MAX_SIZE; i++) vjp_ad += X_dir[i] * X_b[i][idir];
        {
            float abs_err = fabs(vjp_fd - vjp_ad);
            float ref = (fabs(vjp_ad) > 1e-10f) ? fabs(vjp_ad) : 1e-10f;
            float bound = atol + rtol * ref;
            if (abs_err > bound) has_large_errors = 1;
            { float r = abs_err / bound; if (r > max_error) max_error = r; }
        }
    }
    printf("Maximum error ratio: %.6e\n", (double)max_error);
    if (has_large_errors) { printf("FAIL: Large errors in derivatives\n"); return 1; }
    if (max_error < 0.5f) { printf("PASS: Derivatives accurate to machine precision\n"); return 0; }
    if (max_error < 1.0f) { printf("PASS: Derivatives reasonably accurate\n"); return 0; }
    printf("WARNING: Derivatives may have significant errors\n"); return 0;
}

