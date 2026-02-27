/* Test program for cblas_chbmv forward vector (dv) differentiation */
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

extern void cblas_chbmv_dv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, CBLAS_INT N, CBLAS_INT K, const void *alpha, const void *alphad, const void *A, void *Ad, CBLAS_INT lda, const void *X, void *Xd, CBLAS_INT incX, const void *beta, const void *betad, void *Y, void *Yd, CBLAS_INT incY, int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match _d test) */
    float atol = 5.0e-3f, rtol = 5.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| (slightly looser than _d for multi-direction FD) */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements (same as _d) */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO Uplo = CblasUpper;
    CBLAS_INT N = TEST_SIZE;
    CBLAS_INT K = (TEST_SIZE > 1) ? TEST_SIZE - 1 : 0;  /* band width: LDA >= K+1 (match BLAS/test) */
    float complex alpha;
    float complex alphad[NBDirsMax];
    float complex alpha_orig;
    float complex alphad_orig[NBDirsMax];
    float complex A[MAX_SIZE * MAX_SIZE];
    float complex Ad[MAX_SIZE * MAX_SIZE][NBDirsMax];
    float complex A_orig[MAX_SIZE * MAX_SIZE];
    float complex Ad_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT lda = MAX_SIZE;
    float complex X[MAX_SIZE];
    float complex Xd[MAX_SIZE][NBDirsMax];
    float complex X_orig[MAX_SIZE];
    float complex Xd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incX = 1;
    float complex beta;
    float complex betad[NBDirsMax];
    float complex beta_orig;
    float complex betad_orig[NBDirsMax];
    float complex Y[MAX_SIZE];
    float complex Yd[MAX_SIZE][NBDirsMax];
    float complex Y_orig[MAX_SIZE];
    float complex Yd_orig[MAX_SIZE][NBDirsMax];
    CBLAS_INT incY = 1;
    float complex Y_output[MAX_SIZE];
    float complex Y_ad_output[MAX_SIZE];
    float complex Y_forward[MAX_SIZE], Y_backward[MAX_SIZE];

    /* Initialize test data with random numbers (matching _d and Fortran pattern) */
    srand(42);
    alpha = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    /* A: upper band storage (K+1) x N; full a(i,j) at A[K+i-j + j*lda], i = max(0,j-K)..j */
    memset(A, 0, sizeof(A));
    for (j = 0; j < MAX_SIZE; j++) {
        int first_row = (j >= K) ? (j - K) : 0;  /* full row i from first_row..j */
        for (i = first_row; i <= j && i < MAX_SIZE; i++) {
            int band_row = K + i - j;  /* BLAS: a(i,j) -> A(band_row, j) */
            if (i == j) {  /* diagonal: real for Hermitian */
                A[band_row + j * MAX_SIZE] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
            } else {
                A[band_row + j * MAX_SIZE] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f + I * (((float)rand() / RAND_MAX) * 2.0f - 1.0f);
            }
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
    cblas_chbmv(
        layout,
        Uplo,
        N,
        K,
        (const void *)&alpha,
        A,
        lda,
        X,
        incX,
        (const void *)&beta,
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
    cblas_chbmv_dv(
        layout,
        Uplo,
        N,
        K,
        (const void *)&alpha, alphad,
        A, Ad,
        lda,
        X, Xd,
        incX,
        (const void *)&beta, betad,
        Y, Yd,
        incY,
        nbdirs
    );
    memcpy(Y_ad_output, Y, sizeof(Y));

    /* Verify AD primal output matches original (same as _d) */
    {
        float output_diff_max = 0.0f;
        for (i = 0; i < MAX_SIZE; i++) {
            float diff = cabs(Y_ad_output[i] - Y_output[i]);
            if (diff > output_diff_max) output_diff_max = diff;
        }
        if (output_diff_max > 1.0e-10f) {
            printf("WARNING: AD function output differs from original (%s): max_diff=%.6e\n", "Y", (double)output_diff_max);
        }
    }

    /* Compare results using finite differences (same structure as _d) */
    printf("Testing %s differentiation...\n", "cblas_chbmv");
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
        cblas_chbmv(
        layout,
        Uplo,
        N,
        K,
        (const void *)&alpha,
        A,
        lda,
        X,
        incX,
        (const void *)&beta,
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
        cblas_chbmv(
        layout,
        Uplo,
        N,
        K,
        (const void *)&alpha,
        A,
        lda,
        X,
        incX,
        (const void *)&beta,
        Y,
        incY
        );
        memcpy(Y_backward, Y, sizeof(Y));
        /* Central diff vs derivative array(s) */
        for (i = 0; i < MAX_SIZE; i++) {
            float fd = (Y_forward[i] - Y_backward[i]) / (2.0 * h);
            float ad = Yd[i][idir];
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

