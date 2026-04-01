/* Test program for cblas_sgbmv forward vector (dv) differentiation */
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

extern void cblas_sgbmv_dv(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA, CBLAS_INT M, CBLAS_INT N, CBLAS_INT KL, CBLAS_INT KU, float alpha, float alphad[NBDirsMax], const float *A, float (*Ad)[NBDirsMax], CBLAS_INT lda, const float *X, float (*Xd)[NBDirsMax], CBLAS_INT incX, float beta, float betad[NBDirsMax], float *Y, float (*Yd)[NBDirsMax], CBLAS_INT incY, int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match _d test) */
    float atol = 5.0e-3f, rtol = 5.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| (slightly looser than _d for multi-direction FD) */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements (same as _d) */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE TransA = CblasNoTrans;
    CBLAS_INT M = TEST_SIZE;
    CBLAS_INT N = TEST_SIZE;
    CBLAS_INT KL = 1;  /* band width: LDA >= KL+KU+1 (match BLAS/test) */
    CBLAS_INT KU = 1;  /* band width: LDA >= KL+KU+1 (match BLAS/test) */
    float alpha;
    float alphad[NBDirsMax];
    float alpha_orig;
    float alphad_orig[NBDirsMax];
    float A[MAX_SIZE * MAX_SIZE];
    float Ad[MAX_SIZE * MAX_SIZE][NBDirsMax];
    float A_orig[MAX_SIZE * MAX_SIZE];
    float Ad_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT lda = MAX_SIZE;
    float X[MAX_SIZE];
    float Xd[MAX_SIZE][NBDirsMax];
    float X_orig[MAX_SIZE];
    float Xd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incX = 1;
    float beta;
    float betad[NBDirsMax];
    float beta_orig;
    float betad_orig[NBDirsMax];
    float Y[MAX_SIZE];
    float Yd[MAX_SIZE][NBDirsMax];
    float Y_orig[MAX_SIZE];
    float Yd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incY = 1;
    float Y_output[MAX_SIZE];
    float Y_ad_output[MAX_SIZE];
    float Y_forward[MAX_SIZE], Y_backward[MAX_SIZE];

    /* Initialize test data with random numbers (matching _d and Fortran pattern) */
    srand(42);
    alpha = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    /* A: general band storage (KL+KU+1) x N (match BLAS/test) */
    memset(A, 0, sizeof(A));
    for (j = 0; j < MAX_SIZE; j++) {
        int band_rows = KL + KU + 1;
        for (i = 0; i < band_rows; i++) {
            A[i + j * MAX_SIZE] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
        }
    }
    for (i = 0; i < MAX_SIZE; i++) X[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    beta = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) Y[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    /* Initialize derivative seeds (match _d order) */
    for (idir = 0; idir < NBDirsMax; idir++) alphad[idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Ad[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Xd[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (idir = 0; idir < NBDirsMax; idir++) betad[idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Yd[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;

    /* Store originals */
    alpha_orig = alpha;
    memcpy(alphad_orig, alphad, sizeof(alphad));
    memcpy(A_orig, A, sizeof(A));
    memcpy(Ad_orig, Ad, sizeof(Ad));
    memcpy(X_orig, X, sizeof(X));
    memcpy(Xd_orig, Xd, sizeof(Xd));
    beta_orig = beta;
    memcpy(betad_orig, betad, sizeof(betad));
    memcpy(Y_orig, Y, sizeof(Y));
    memcpy(Yd_orig, Yd, sizeof(Yd));

    /* Warmup + primal call, save output(s) */
    cblas_sgbmv(
        layout,
        TransA,
        M,
        N,
        KL,
        KU,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY
    );
    memcpy(Y_output, Y, sizeof(Y));

    /* Restore all inputs and derivative seeds */
    alpha = alpha_orig;
    memcpy(alphad, alphad_orig, sizeof(alphad));
    memcpy(A, A_orig, sizeof(A));
    memcpy(Ad, Ad_orig, sizeof(Ad));
    memcpy(X, X_orig, sizeof(X));
    memcpy(Xd, Xd_orig, sizeof(Xd));
    beta = beta_orig;
    memcpy(betad, betad_orig, sizeof(betad));
    memcpy(Y, Y_orig, sizeof(Y));
    memcpy(Yd, Yd_orig, sizeof(Yd));

    /* Call _dv (implementation uses const void* for alpha/beta in complex, so pass pointers) */
    cblas_sgbmv_dv(
        layout,
        TransA,
        M,
        N,
        KL,
        KU,
        alpha, alphad,
        A, Ad,
        lda,
        X, Xd,
        incX,
        beta, betad,
        Y, Yd,
        incY,
        nbdirs
    );
    memcpy(Y_ad_output, Y, sizeof(Y));

    /* Verify AD primal output matches original (same as _d) */
    {
        float output_diff_max = 0.0f;
        for (i = 0; i < MAX_SIZE; i++) {
            float diff = fabs(Y_ad_output[i] - Y_output[i]);
            if (diff > output_diff_max) output_diff_max = diff;
        }
        if (output_diff_max > 1.0e-10f) {
            printf("WARNING: AD function output differs from original (%s): max_diff=%.6e\n", "Y", (double)output_diff_max);
        }
    }

    /* Compare results using finite differences (same structure as _d) */
    printf("Testing %s differentiation...\n", "cblas_sgbmv");
    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals (matching _d) */
        alpha = alpha_orig;
        memcpy(A, A_orig, sizeof(A));
        memcpy(X, X_orig, sizeof(X));
        beta = beta_orig;
        memcpy(Y, Y_orig, sizeof(Y));
        /* Forward perturbation: x + h * x_d (same order as _d) */
        alpha += h * alphad_orig[idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) A[j] += h * Ad_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) X[j] += h * Xd_orig[j][idir];
        beta += h * betad_orig[idir];
        for (j = 0; j < MAX_SIZE; j++) Y[j] += h * Yd_orig[j][idir];
        cblas_sgbmv(
        layout,
        TransA,
        M,
        N,
        KL,
        KU,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY
        );
        memcpy(Y_forward, Y, sizeof(Y));
        /* Restore primals (matching _d) */
        alpha = alpha_orig;
        memcpy(A, A_orig, sizeof(A));
        memcpy(X, X_orig, sizeof(X));
        beta = beta_orig;
        memcpy(Y, Y_orig, sizeof(Y));
        /* Backward perturbation: x - h * x_d (same order as _d) */
        alpha -= h * alphad_orig[idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) A[j] -= h * Ad_orig[j][idir];
        for (j = 0; j < MAX_SIZE; j++) X[j] -= h * Xd_orig[j][idir];
        beta -= h * betad_orig[idir];
        for (j = 0; j < MAX_SIZE; j++) Y[j] -= h * Yd_orig[j][idir];
        cblas_sgbmv(
        layout,
        TransA,
        M,
        N,
        KL,
        KU,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY
        );
        memcpy(Y_backward, Y, sizeof(Y));
        /* Central diff vs derivative array(s) */
        for (i = 0; i < MAX_SIZE; i++) {
            float fd = (Y_forward[i] - Y_backward[i]) / (2.0 * h);
            float ad = Yd[i][idir];
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

