/* Test program for cblas_sscal forward vector (dv) differentiation */
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
extern void set_isize1ofsx_(int *val);

extern void cblas_sscal_dv(CBLAS_INT N, float alpha, float alphad[NBDirsMax], float *X, float (*Xd)[NBDirsMax], CBLAS_INT incX, int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofsx_(&diffblas_isize);
    }
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match _d test) */
    float atol = 5.0e-3f, rtol = 5.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| (slightly looser than _d for multi-direction FD) */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements (same as _d) */

    CBLAS_INT N = TEST_SIZE;
    float alpha;
    float alphad[NBDirsMax];
    float alpha_orig;
    float alphad_orig[NBDirsMax];
    float X[MAX_SIZE];
    float Xd[MAX_SIZE][NBDirsMax];
    float X_orig[MAX_SIZE];
    float Xd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incX = 1;
    float X_output[MAX_SIZE];
    float X_ad_output[MAX_SIZE];
    float X_forward[MAX_SIZE], X_backward[MAX_SIZE];

    /* Initialize test data with random numbers (matching _d and Fortran pattern) */
    srand(42);
    alpha = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) X[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    /* Initialize derivative seeds (match _d order) */
    for (idir = 0; idir < NBDirsMax; idir++) alphad[idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Xd[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;

    /* Store originals */
    alpha_orig = alpha;
    memcpy(alphad_orig, alphad, sizeof(alphad));
    memcpy(X_orig, X, sizeof(X));
    memcpy(Xd_orig, Xd, sizeof(Xd));

    /* Warmup + primal call, save output(s) */
    cblas_sscal(
        N,
        alpha,
        X,
        incX
    );
    memcpy(X_output, X, sizeof(X));

    /* Restore all inputs and derivative seeds */
    alpha = alpha_orig;
    memcpy(alphad, alphad_orig, sizeof(alphad));
    memcpy(X, X_orig, sizeof(X));
    memcpy(Xd, Xd_orig, sizeof(Xd));

    /* Call _dv (implementation uses const void* for alpha/beta in complex, so pass pointers) */
    cblas_sscal_dv(
        N,
        alpha, alphad,
        X, Xd,
        incX,
        nbdirs
    );
    memcpy(X_ad_output, X, sizeof(X));

    /* Verify AD primal output matches original (same as _d) */
    {
        float output_diff_max = 0.0f;
        for (i = 0; i < MAX_SIZE; i++) {
            float diff = fabs(X_ad_output[i] - X_output[i]);
            if (diff > output_diff_max) output_diff_max = diff;
        }
        if (output_diff_max > 1.0e-10f) {
            printf("WARNING: AD function output differs from original (%s): max_diff=%.6e\n", "X", (double)output_diff_max);
        }
    }

    /* Compare results using finite differences (same structure as _d) */
    printf("Testing %s differentiation...\n", "cblas_sscal");
    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals (matching _d) */
        alpha = alpha_orig;
        memcpy(X, X_orig, sizeof(X));
        /* Forward perturbation: x + h * x_d (same order as _d) */
        alpha += h * alphad_orig[idir];
        for (j = 0; j < MAX_SIZE; j++) X[j] += h * Xd_orig[j][idir];
        cblas_sscal(
        N,
        alpha,
        X,
        incX
        );
        memcpy(X_forward, X, sizeof(X));
        /* Restore primals (matching _d) */
        alpha = alpha_orig;
        memcpy(X, X_orig, sizeof(X));
        /* Backward perturbation: x - h * x_d (same order as _d) */
        alpha -= h * alphad_orig[idir];
        for (j = 0; j < MAX_SIZE; j++) X[j] -= h * Xd_orig[j][idir];
        cblas_sscal(
        N,
        alpha,
        X,
        incX
        );
        memcpy(X_backward, X, sizeof(X));
        /* Central diff vs derivative array(s) */
        for (i = 0; i < MAX_SIZE; i++) {
            float fd = (X_forward[i] - X_backward[i]) / (2.0 * h);
            float ad = Xd[i][idir];
            float abs_err = fabs(fd - ad);
            float ad_ref = (fabs(ad) > 1e-10) ? fabs(ad) : 1e-10;
            float bound = atol + rtol * ad_ref;
            if (abs_err > bound) { has_large_errors = 1; }
            float r = abs_err / bound;
            if (r > max_error) max_error = r;
        }
    }
    printf("Maximum error ratio (abs_error/error_bound): %.6e\n", (double)max_error);
    if (has_large_errors) {
        printf("FAIL: Large errors detected in derivatives\n");
        return 1;
    }
    else if (max_error < 0.5f) {
        printf("PASS: Derivatives are accurate to machine precision\n");
        return 0;
    }
    else if (max_error < 2.0f) {
        printf("PASS: Derivatives are reasonably accurate\n");
        return 0;
    } else {
        printf("WARNING: Derivatives may have significant errors\n");
        return 0;
    }
}

