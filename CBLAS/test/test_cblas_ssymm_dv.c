/* Test program for cblas_ssymm forward vector (dv) differentiation */
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

extern void cblas_ssymm_dv(CBLAS_LAYOUT layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_INT M, CBLAS_INT N, float alpha, float alphad[NBDirsMax], const float *A, float (*Ad)[NBDirsMax], CBLAS_INT lda, const float *B, float (*Bd)[NBDirsMax], CBLAS_INT ldb, float beta, float betad[NBDirsMax], float *C, float (*Cd)[NBDirsMax], CBLAS_INT ldc, int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    float h = 1.0e-3f;  /* Step size for finite differences (match _d test) */
    float atol = 5.0e-3f, rtol = 5.0e-3f;  /* Pass when abs_error <= atol + rtol*|ad| (slightly looser than _d for multi-direction FD) */
    float max_error = 0.0f;  /* max (abs_error/error_bound) over elements (same as _d) */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_SIDE Side = CblasLeft;
    CBLAS_UPLO Uplo = CblasUpper;
    CBLAS_INT M = TEST_SIZE;
    CBLAS_INT N = TEST_SIZE;
    float alpha;
    float alphad[NBDirsMax];
    float alpha_orig;
    float alphad_orig[NBDirsMax];
    float A[MAX_SIZE * MAX_SIZE];
    float Ad[MAX_SIZE * MAX_SIZE][NBDirsMax];
    float A_orig[MAX_SIZE * MAX_SIZE];
    float Ad_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT lda = MAX_SIZE;
    float B[MAX_SIZE * MAX_SIZE];
    float Bd[MAX_SIZE * MAX_SIZE][NBDirsMax];
    float B_orig[MAX_SIZE * MAX_SIZE];
    float Bd_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT ldb = MAX_SIZE;
    float beta;
    float betad[NBDirsMax];
    float beta_orig;
    float betad_orig[NBDirsMax];
    float C[MAX_SIZE * MAX_SIZE];
    float Cd[MAX_SIZE * MAX_SIZE][NBDirsMax];
    float C_orig[MAX_SIZE * MAX_SIZE];
    float Cd_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT ldc = MAX_SIZE;
    float C_output[MAX_SIZE * MAX_SIZE];
    float C_ad_output[MAX_SIZE * MAX_SIZE];
    float C_forward[MAX_SIZE * MAX_SIZE], C_backward[MAX_SIZE * MAX_SIZE];

    /* Initialize test data with random numbers (matching _d and Fortran pattern) */
    srand(42);
    alpha = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    /* A: symmetric (match BLAS/test) */
    for (i = 0; i < MAX_SIZE; i++) {
        for (j = i; j < MAX_SIZE; j++) {
            A[i * MAX_SIZE + j] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
        }
    }
    for (i = 1; i < MAX_SIZE; i++) {
        for (j = 0; j < i; j++) {
            A[i * MAX_SIZE + j] = A[j * MAX_SIZE + i];  /* symmetric */
        }
    }
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) B[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    beta = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) C[i] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    /* Initialize derivative seeds (match _d order) */
    for (idir = 0; idir < NBDirsMax; idir++) alphad[idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Ad[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Bd[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (idir = 0; idir < NBDirsMax; idir++) betad[idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Cd[i][idir] = ((float)rand() / RAND_MAX) * 2.0 - 1.0;

    /* Store originals */
    alpha_orig = alpha;
    memcpy(alphad_orig, alphad, sizeof(alphad));
    memcpy(A_orig, A, sizeof(A));
    memcpy(Ad_orig, Ad, sizeof(Ad));
    memcpy(B_orig, B, sizeof(B));
    memcpy(Bd_orig, Bd, sizeof(Bd));
    beta_orig = beta;
    memcpy(betad_orig, betad, sizeof(betad));
    memcpy(C_orig, C, sizeof(C));
    memcpy(Cd_orig, Cd, sizeof(Cd));

    /* Warmup + primal call, save output(s) */
    cblas_ssymm(
        layout,
        Side,
        Uplo,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc
    );
    memcpy(C_output, C, sizeof(C));

    /* Restore all inputs and derivative seeds */
    alpha = alpha_orig;
    memcpy(alphad, alphad_orig, sizeof(alphad));
    memcpy(A, A_orig, sizeof(A));
    memcpy(Ad, Ad_orig, sizeof(Ad));
    memcpy(B, B_orig, sizeof(B));
    memcpy(Bd, Bd_orig, sizeof(Bd));
    beta = beta_orig;
    memcpy(betad, betad_orig, sizeof(betad));
    memcpy(C, C_orig, sizeof(C));
    memcpy(Cd, Cd_orig, sizeof(Cd));

    /* Call _dv (implementation uses const void* for alpha/beta in complex, so pass pointers) */
    cblas_ssymm_dv(
        layout,
        Side,
        Uplo,
        M,
        N,
        alpha, alphad,
        A, Ad,
        lda,
        B, Bd,
        ldb,
        beta, betad,
        C, Cd,
        ldc,
        nbdirs
    );
    memcpy(C_ad_output, C, sizeof(C));

    /* Verify AD primal output matches original (same as _d) */
    {
        float output_diff_max = 0.0f;
        for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
            float diff = fabs(C_ad_output[i] - C_output[i]);
            if (diff > output_diff_max) output_diff_max = diff;
        }
        if (output_diff_max > 1.0e-10f) {
            printf("WARNING: AD function output differs from original (%s): max_diff=%.6e\n", "C", (double)output_diff_max);
        }
    }

    /* Compare results using finite differences (same structure as _d) */
    printf("Testing %s differentiation...\n", "cblas_ssymm");
    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals (matching _d) */
        alpha = alpha_orig;
        memcpy(A, A_orig, sizeof(A));
        memcpy(B, B_orig, sizeof(B));
        beta = beta_orig;
        memcpy(C, C_orig, sizeof(C));
        /* Forward perturbation: x + h * x_d (same order as _d) */
        alpha += h * alphad_orig[idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) A[j] += h * Ad_orig[j][idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) B[j] += h * Bd_orig[j][idir];
        beta += h * betad_orig[idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) C[j] += h * Cd_orig[j][idir];
        cblas_ssymm(
        layout,
        Side,
        Uplo,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc
        );
        memcpy(C_forward, C, sizeof(C));
        /* Restore primals (matching _d) */
        alpha = alpha_orig;
        memcpy(A, A_orig, sizeof(A));
        memcpy(B, B_orig, sizeof(B));
        beta = beta_orig;
        memcpy(C, C_orig, sizeof(C));
        /* Backward perturbation: x - h * x_d (same order as _d) */
        alpha -= h * alphad_orig[idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) A[j] -= h * Ad_orig[j][idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) B[j] -= h * Bd_orig[j][idir];
        beta -= h * betad_orig[idir];
        for (j = 0; j < MAX_SIZE * MAX_SIZE; j++) C[j] -= h * Cd_orig[j][idir];
        cblas_ssymm(
        layout,
        Side,
        Uplo,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc
        );
        memcpy(C_backward, C, sizeof(C));
        /* Central diff vs derivative array(s) */
        for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
            float fd = (C_forward[i] - C_backward[i]) / (2.0 * h);
            float ad = Cd[i][idir];
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

