/* Test program for cblas_dsyr2k vector reverse mode (VJP verification, generic, loop over directions) */
/* Generated automatically by run_tapenade_cblas.py */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cblas.h"
#include "cblas_f77.h"
#include "cblas_bv.h"

#ifndef NBDirsMax
#define NBDirsMax 4
#endif
#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE

static int compare_abs_d(const void *a, const void *b) { double x = fabs(*(const double*)a), y = fabs(*(const double*)b); return (x > y) - (x < y); }

extern void cblas_dsyr2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE Trans, const CBLAS_INT N, const CBLAS_INT K, const double alpha, const double *A, const CBLAS_INT lda, const double *B, const CBLAS_INT ldb, const double beta, double *C, const CBLAS_INT ldc);
/* cblas_*_bv from cblas_bv.h */

int main(void) {
    int i, j, idx, idir, nbdirs = NBDirsMax, n_products;
    int has_large_errors = 0;
    double h = 1.0e-7;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double max_error = 0.0;
    double vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT k = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT ldb = MAX_SIZE;
    CBLAS_INT ldc = MAX_SIZE;

    double alpha, alpha_b[NBDirsMax], alpha_orig, alpha_dir, alpha_b_orig[NBDirsMax];
    double A[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    double A_b[MAX_SIZE*MAX_SIZE][NBDirsMax], A_b_orig[MAX_SIZE*MAX_SIZE][NBDirsMax];
    double B[MAX_SIZE*MAX_SIZE], B_orig[MAX_SIZE*MAX_SIZE], B_dir[MAX_SIZE*MAX_SIZE];
    double B_b[MAX_SIZE*MAX_SIZE][NBDirsMax], B_b_orig[MAX_SIZE*MAX_SIZE][NBDirsMax];
    double beta, beta_b[NBDirsMax], beta_orig, beta_dir, beta_b_orig[NBDirsMax];
    double C[MAX_SIZE*MAX_SIZE], C_orig[MAX_SIZE*MAX_SIZE], C_dir[MAX_SIZE*MAX_SIZE];
    double C_b[MAX_SIZE*MAX_SIZE][NBDirsMax], C_b_orig[MAX_SIZE*MAX_SIZE][NBDirsMax];
    double C_plus[MAX_SIZE*MAX_SIZE], C_minus[MAX_SIZE*MAX_SIZE];

    srand(42);
    alpha = (rand()/(double)RAND_MAX)*2.0 - 1.0; alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { B[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { B_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    beta = (rand()/(double)RAND_MAX)*2.0 - 1.0; beta_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { C[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { C_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }

    alpha_orig = alpha;
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
    memcpy(B_orig, B, sizeof(B[0])*(MAX_SIZE*MAX_SIZE));
    beta_orig = beta;
    memcpy(C_orig, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));

    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) { C_b[i][j] = (rand()/(double)RAND_MAX)*2.0 - 1.0; C_b_orig[i][j] = C_b[i][j]; }
    for (j = 0; j < NBDirsMax; j++) alpha_b[j] = 0.0;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) A_b[i][j] = 0.0;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) B_b[i][j] = 0.0;
    for (j = 0; j < NBDirsMax; j++) beta_b[j] = 0.0;

    cblas_dsyr2k_bv(layout, uplo, trans, n, k, alpha, &alpha_b, A, A_b, lda, B, &B_b[0][0], ldb, beta, &beta_b, C, &C_b[0][0], ldc, nbdirs);

    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals for this direction */
        alpha = alpha_orig;
        memcpy(A, A_orig, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
        memcpy(B, B_orig, sizeof(B[0])*(MAX_SIZE*MAX_SIZE));
        beta = beta_orig;
        memcpy(C, C_orig, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));
        /* Random direction for this idir */
        alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) B_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        beta_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) C_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        /* Forward */
        alpha = alpha_orig + h * alpha_dir;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] + h * A_dir[i];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) B[i] = B_orig[i] + h * B_dir[i];
        beta = beta_orig + h * beta_dir;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) C[i] = C_orig[i] + h * C_dir[i];
        cblas_dsyr2k(layout, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        memcpy(C_plus, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));
        /* Backward */
        alpha = alpha_orig - h * alpha_dir;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) B[i] = B_orig[i] - h * B_dir[i];
        beta = beta_orig - h * beta_dir;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) C[i] = C_orig[i] - h * C_dir[i];
        cblas_dsyr2k(layout, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        memcpy(C_minus, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));

        vjp_fd = 0.0;
        {
            double temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = C_b_orig[i][idir] * ((C_plus[i] - C_minus[i]) / (2.0*h));
            qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(double), compare_abs_d);
            for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_fd += temp_products[idx];
        }
        vjp_ad = 0.0;
        vjp_ad += alpha_dir * alpha_b[idir];
        {
            double temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = A_dir[i] * A_b[i][idir];
            qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(double), compare_abs_d);
            for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }
        {
            double temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = B_dir[i] * B_b[i][idir];
            qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(double), compare_abs_d);
            for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }
        vjp_ad += beta_dir * beta_b[idir];
        {
            double temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = C_dir[i] * C_b[i][idir];
            qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(double), compare_abs_d);
            for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }

        {
            double abs_err = fabs(vjp_fd - vjp_ad);
            double abs_reference = fabs(vjp_ad);
            double error_bound = atol + rtol * (abs_reference > 1e-10 ? abs_reference : 1e-10);
            if (abs_err > error_bound) has_large_errors = 1;
            { double r = abs_err / error_bound; if (r > max_error) max_error = r; }
        }
    }

    printf("Maximum error ratio: %.6e\n", (double)max_error);
    if (has_large_errors) { printf("FAIL: Large errors in derivatives\n"); return 1; }
    if (max_error < 0.5) { printf("PASS: Derivatives accurate to machine precision\n"); return 0; }
    if (max_error < 1.0) { printf("PASS: Derivatives reasonably accurate\n"); return 0; }
    printf("WARNING: Derivatives may have significant errors\n"); return 0;
}

