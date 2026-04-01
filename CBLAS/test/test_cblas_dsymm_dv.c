/* Test program for cblas_dsymm forward vector (dv) differentiation */
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

extern void cblas_dsymm_dv(CBLAS_LAYOUT layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_INT M, CBLAS_INT N, double alpha, double alphad[NBDirsMax], const double *A, double (*Ad)[NBDirsMax], CBLAS_INT lda, const double *B, double (*Bd)[NBDirsMax], CBLAS_INT ldb, double beta, double betad[NBDirsMax], double *C, double (*Cd)[NBDirsMax], CBLAS_INT ldc, int nbdirs);

int main(void) {
    int i, j, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    double h = 1.0e-6;  /* Step size for finite differences (match _d test) */
    double atol = 1.0e-5, rtol = 1.0e-5;  /* Pass when abs_error <= atol + rtol*|ad| (same as _d) */
    double max_error = 0.0;  /* max (abs_error/error_bound) over elements (same as _d) */

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_SIDE Side = CblasLeft;
    CBLAS_UPLO Uplo = CblasUpper;
    CBLAS_INT M = TEST_SIZE;
    CBLAS_INT N = TEST_SIZE;
    double alpha;
    double alphad[NBDirsMax];
    double alpha_orig;
    double alphad_orig[NBDirsMax];
    double A[MAX_SIZE * MAX_SIZE];
    double Ad[MAX_SIZE * MAX_SIZE][NBDirsMax];
    double A_orig[MAX_SIZE * MAX_SIZE];
    double Ad_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT lda = MAX_SIZE;
    double B[MAX_SIZE * MAX_SIZE];
    double Bd[MAX_SIZE * MAX_SIZE][NBDirsMax];
    double B_orig[MAX_SIZE * MAX_SIZE];
    double Bd_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT ldb = MAX_SIZE;
    double beta;
    double betad[NBDirsMax];
    double beta_orig;
    double betad_orig[NBDirsMax];
    double C[MAX_SIZE * MAX_SIZE];
    double Cd[MAX_SIZE * MAX_SIZE][NBDirsMax];
    double C_orig[MAX_SIZE * MAX_SIZE];
    double Cd_orig[MAX_SIZE * MAX_SIZE][NBDirsMax];
    CBLAS_INT ldc = MAX_SIZE;
    double C_output[MAX_SIZE * MAX_SIZE];
    double C_ad_output[MAX_SIZE * MAX_SIZE];
    double C_forward[MAX_SIZE * MAX_SIZE], C_backward[MAX_SIZE * MAX_SIZE];

    /* Initialize test data with random numbers (matching _d and Fortran pattern) */
    srand(42);
    alpha = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    /* A: symmetric (match BLAS/test) */
    for (i = 0; i < MAX_SIZE; i++) {
        for (j = i; j < MAX_SIZE; j++) {
            A[i * MAX_SIZE + j] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        }
    }
    for (i = 1; i < MAX_SIZE; i++) {
        for (j = 0; j < i; j++) {
            A[i * MAX_SIZE + j] = A[j * MAX_SIZE + i];  /* symmetric */
        }
    }
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) B[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    beta = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) C[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    /* Initialize derivative seeds (match _d order) */
    for (idir = 0; idir < NBDirsMax; idir++) alphad[idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Ad[i][idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Bd[i][idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (idir = 0; idir < NBDirsMax; idir++) betad[idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) for (idir = 0; idir < NBDirsMax; idir++) Cd[i][idir] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;

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
    cblas_dsymm(
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
    cblas_dsymm_dv(
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
        double output_diff_max = 0.0;
        for (i = 0; i < MAX_SIZE * MAX_SIZE; i++) {
            double diff = fabs(C_ad_output[i] - C_output[i]);
            if (diff > output_diff_max) output_diff_max = diff;
        }
        if (output_diff_max > 1.0e-10) {
            printf("WARNING: AD function output differs from original (%s): max_diff=%.6e\n", "C", (double)output_diff_max);
        }
    }

    /* Compare results using finite differences (same structure as _d) */
    printf("Testing %s differentiation...\n", "cblas_dsymm");
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
        cblas_dsymm(
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
        cblas_dsymm(
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
            double fd = (C_forward[i] - C_backward[i]) / (2.0 * h);
            double ad = Cd[i][idir];
            double abs_err = fabs(fd - ad);
            double ad_ref = (fabs(ad) > 1e-10) ? fabs(ad) : 1e-10;
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
    else if (max_error < 0.5) {
        printf("PASS: Derivatives are accurate to machine precision\n");
        return 0;
    }
    else if (max_error < 1.0) {
        printf("PASS: Derivatives are reasonably accurate\n");
        return 0;
    } else {
        printf("WARNING: Derivatives may have significant errors\n");
        return 0;
    }
}

