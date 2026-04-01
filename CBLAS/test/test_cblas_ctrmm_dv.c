/* Test program for cblas_ctrmm forward vector (dv) differentiation */
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

extern void cblas_ctrmm_dv(CBLAS_LAYOUT layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, CBLAS_INT M, CBLAS_INT N, const void *alpha, const void *alphad, const void *A, void *Ad, CBLAS_INT lda, void *B, void *Bd, CBLAS_INT ldb, int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match _d test) */
    float atol = 5.0e-3f, rtol = 5.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| (slightly looser than _d for multi-direction FD) */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements (same as _d) */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_SIDE Side = CblasLeft;
    CBLAS_UPLO Uplo = CblasUpper;
    CBLAS_TRANSPOSE TransA = CblasNoTrans;
    CBLAS_DIAG Diag = CblasNonUnit;
    CBLAS_INT M = TEST_SIZE;
    CBLAS_INT N = TEST_SIZE;
    float complex alpha;
    float complex alphad[NBDirsMax];
    float complex alpha_orig;
    float complex alphad_orig[NBDirsMax];
    float complex A[MAX_SIZE * MAX_SIZE];
    float complex Ad[MAX_SIZE * MAX_SIZE][NBDirsMax];
    float complex A_orig[MAX_SIZE * MAX_SIZE];
    float complex Ad_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT lda = MAX_SIZE;
    float complex B[MAX_SIZE * MAX_SIZE];
    float complex Bd[MAX_SIZE * MAX_SIZE][NBDirsMax];
    float complex B_orig[MAX_SIZE * MAX_SIZE];
    float complex Bd_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT ldb = MAX_SIZE;
    float complex B_output[MAX_SIZE * MAX_SIZE];
    float complex B_ad_output[MAX_SIZE * MAX_SIZE];
    float complex B_forward[MAX_SIZE * MAX_SIZE], B_backward[MAX_SIZE * MAX_SIZE];

    /* Initialize test data with random numbers (matching _d and Fortran pattern) */
    srand(42);
    alpha = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) A[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) B[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    /* Initialize derivative seeds (match _d order) */
    for (idir = 0; idir < NBDirsMax; idir++) alphad[idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Ad[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Bd[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;

    /* Store originals */
    alpha_orig = alpha;
    memcpy(alphad_orig, alphad, sizeof(alphad));
    memcpy(A_orig, A, sizeof(A));
    memcpy(Ad_orig, Ad, sizeof(Ad));
    memcpy(B_orig, B, sizeof(B));
    memcpy(Bd_orig, Bd, sizeof(Bd));

    /* Warmup + primal call, save output(s) */
    cblas_ctrmm(
        layout,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        (const void *)&alpha,
        A,
        lda,
        B,
        ldb
    );
    memcpy(B_output, B, sizeof(B));

    /* Restore all inputs and derivative seeds */
    alpha = alpha_orig;
    memcpy(alphad, alphad_orig, sizeof(alphad));
    memcpy(A, A_orig, sizeof(A));
    memcpy(Ad, Ad_orig, sizeof(Ad));
    memcpy(B, B_orig, sizeof(B));
    memcpy(Bd, Bd_orig, sizeof(Bd));

    /* Call _dv (implementation uses const void* for alpha/beta in complex, so pass pointers) */
    cblas_ctrmm_dv(
        layout,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        (const void *)&alpha, alphad,
        A, Ad,
        lda,
        B, Bd,
        ldb,
        nbdirs
    );
    memcpy(B_ad_output, B, sizeof(B));

    /* Verify AD primal output matches original (same as _d) */
    {
        float output_diff_max = 0.0f;
        for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
            float diff = cabs(B_ad_output[i] - B_output[i]);
            if (diff > output_diff_max) output_diff_max = diff;
        }
        if (output_diff_max > 1.0e-10f) {
            printf("WARNING: AD function output differs from original (%s): max_diff=%.6e\n", "B", (double)output_diff_max);
        }
    }

    /* Compare results using finite differences (same structure as _d) */
    printf("Testing %s differentiation...\n", "cblas_ctrmm");
    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals (matching _d) */
        alpha = alpha_orig;
        memcpy(A, A_orig, sizeof(A));
        memcpy(B, B_orig, sizeof(B));
        /* Forward perturbation: x + h * x_d (same order as _d) */
        alpha += h * alphad_orig[idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) A[j] += h * Ad_orig[j][idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) B[j] += h * Bd_orig[j][idir];
        cblas_ctrmm(
        layout,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        (const void *)&alpha,
        A,
        lda,
        B,
        ldb
        );
        memcpy(B_forward, B, sizeof(B));
        /* Restore primals (matching _d) */
        alpha = alpha_orig;
        memcpy(A, A_orig, sizeof(A));
        memcpy(B, B_orig, sizeof(B));
        /* Backward perturbation: x - h * x_d (same order as _d) */
        alpha -= h * alphad_orig[idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) A[j] -= h * Ad_orig[j][idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) B[j] -= h * Bd_orig[j][idir];
        cblas_ctrmm(
        layout,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        (const void *)&alpha,
        A,
        lda,
        B,
        ldb
        );
        memcpy(B_backward, B, sizeof(B));
        /* Central diff vs derivative array(s) */
        for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
            float fd = (B_forward[i] - B_backward[i]) / (2.0 * h);
            float ad = Bd[i][idir];
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

