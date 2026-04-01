/* Test program for cblas_sspr forward vector (dv) differentiation */
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
extern void set_isize1ofap_(int *val);

extern void cblas_sspr_dv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, CBLAS_INT N, float alpha, float alphad[NBDirsMax], const float *X, float (*Xd)[NBDirsMax], CBLAS_INT incX, float *Ap, float (*Apd)[NBDirsMax], int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofap_(&diffblas_isize);
    }
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match _d test) */
    float atol = 5.0e-3f, rtol = 5.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| (slightly looser than _d for multi-direction FD) */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements (same as _d) */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO Uplo = CblasUpper;
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
    float Ap[PACKED_SIZE];
    float Apd[PACKED_SIZE][NBDirsMax];
    float Ap_orig[PACKED_SIZE];
    float Apd_orig[PACKED_SIZE][NBDirsMax];
    float Ap_output[PACKED_SIZE];
    float Ap_ad_output[PACKED_SIZE];
    float Ap_forward[PACKED_SIZE], Ap_backward[PACKED_SIZE];

    /* Initialize test data with random numbers (matching _d and Fortran pattern) */
    srand(42);
    alpha = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) X[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (j = 0; j < MAX_SIZE; j++)
        for (i = 0; i <= j; i++)
            Ap[j * (j + 1) / 2 + i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0f;
    /* Initialize derivative seeds (match _d order) */
    for (idir = 0; idir < NBDirsMax; idir++) alphad[idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Xd[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < PACKED_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Apd[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;

    /* Store originals */
    alpha_orig = alpha;
    memcpy(alphad_orig, alphad, sizeof(alphad));
    memcpy(X_orig, X, sizeof(X));
    memcpy(Xd_orig, Xd, sizeof(Xd));
    memcpy(Ap_orig, Ap, sizeof(Ap));
    memcpy(Apd_orig, Apd, sizeof(Apd));

    /* Warmup + primal call, save output(s) */
    cblas_sspr(
        layout,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Ap
    );
    memcpy(Ap_output, Ap, sizeof(Ap));

    /* Restore all inputs and derivative seeds */
    alpha = alpha_orig;
    memcpy(alphad, alphad_orig, sizeof(alphad));
    memcpy(X, X_orig, sizeof(X));
    memcpy(Xd, Xd_orig, sizeof(Xd));
    memcpy(Ap, Ap_orig, sizeof(Ap));
    memcpy(Apd, Apd_orig, sizeof(Apd));

    /* Call _dv (implementation uses const void* for alpha/beta in complex, so pass pointers) */
    cblas_sspr_dv(
        layout,
        Uplo,
        N,
        alpha, alphad,
        X, Xd,
        incX,
        Ap, Apd,
        nbdirs
    );
    memcpy(Ap_ad_output, Ap, sizeof(Ap));

    /* Verify AD primal output matches original (same as _d) */
    {
        float output_diff_max = 0.0f;
        for (i = 0; i < PACKED_SIZE; i++) {
            float diff = fabs(Ap_ad_output[i] - Ap_output[i]);
            if (diff > output_diff_max) output_diff_max = diff;
        }
        if (output_diff_max > 1.0e-10f) {
            printf("WARNING: AD function output differs from original (%s): max_diff=%.6e\n", "Ap", (double)output_diff_max);
        }
    }

    /* Compare results using finite differences (same structure as _d) */
    printf("Testing %s differentiation...\n", "cblas_sspr");
    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals (matching _d) */
        alpha = alpha_orig;
        memcpy(X, X_orig, sizeof(X));
        memcpy(Ap, Ap_orig, sizeof(Ap));
        /* Forward perturbation: x + h * x_d (same order as _d) */
        alpha += h * alphad_orig[idir];
        for (j = 0; j < MAX_SIZE; j++) X[j] += h * Xd_orig[j][idir];
        for (j = 0; j < PACKED_SIZE; j++) Ap[j] += h * Apd_orig[j][idir];
        cblas_sspr(
        layout,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Ap
        );
        memcpy(Ap_forward, Ap, sizeof(Ap));
        /* Restore primals (matching _d) */
        alpha = alpha_orig;
        memcpy(X, X_orig, sizeof(X));
        memcpy(Ap, Ap_orig, sizeof(Ap));
        /* Backward perturbation: x - h * x_d (same order as _d) */
        alpha -= h * alphad_orig[idir];
        for (j = 0; j < MAX_SIZE; j++) X[j] -= h * Xd_orig[j][idir];
        for (j = 0; j < PACKED_SIZE; j++) Ap[j] -= h * Apd_orig[j][idir];
        cblas_sspr(
        layout,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Ap
        );
        memcpy(Ap_backward, Ap, sizeof(Ap));
        /* Central diff vs derivative array(s) */
        for (i = 0; i < PACKED_SIZE; i++) {
            float fd = (Ap_forward[i] - Ap_backward[i]) / (2.0 * h);
            float ad = Apd[i][idir];
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

