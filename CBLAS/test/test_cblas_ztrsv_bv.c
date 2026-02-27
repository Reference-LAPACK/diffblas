/* Test program for cblas_ztrsv vector reverse mode (VJP verification, generic, loop over directions) */
/* Generated automatically by run_tapenade_cblas.py */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "cblas.h"
#include "cblas_f77.h"
#include "cblas_bv.h"

#ifndef NBDirsMax
#define NBDirsMax 4
#endif
#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE

static int compare_abs_d(const void *a, const void *b) { double x = fabs(*(const double*)a), y = fabs(*(const double*)b); return (x > y) - (x < y); }

extern void cblas_ztrsv(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const CBLAS_INT N, const void *A, const CBLAS_INT lda, void *X, const CBLAS_INT incX);
/* cblas_*_bv from cblas_bv.h */

int main(void) {
    int i, j, idx, idir, nbdirs = NBDirsMax, n_products;
    int has_large_errors = 0;
    double h = 1.0e-7;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double max_error = 0.0;
    double vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE transa = CblasNoTrans;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_DIAG diag = CblasNonUnit;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT incX = 1;

    double complex A[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    double complex A_b[MAX_SIZE*MAX_SIZE][NBDirsMax], A_b_orig[MAX_SIZE*MAX_SIZE][NBDirsMax];
    double complex X[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    double complex X_b[MAX_SIZE][NBDirsMax], X_b_orig[MAX_SIZE][NBDirsMax];
    double complex X_plus[MAX_SIZE], X_minus[MAX_SIZE];

    srand(42);
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }

    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));

    /* Triangular A: zero unused triangle and set unit diagonal if needed */
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
            if (uplo == CblasUpper && i > j) { A[i + j*lda] = 0.0; A_dir[i + j*lda] = 0.0; }
            if (uplo == CblasLower && i < j) { A[i + j*lda] = 0.0; A_dir[i + j*lda] = 0.0; }
        }
        if (diag == CblasUnit) { A[j + j*lda] = 1.0; A_dir[j + j*lda] = 0.0; }
    }
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) { X_b[i][j] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); X_b_orig[i][j] = X_b[i][j]; }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) A_b[i][j] = 0.0;

    cblas_ztrsv_bv(layout, uplo, transa, diag, n, A, (void*)A_b, lda, X, (void*)X_b, incX, nbdirs);

    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals for this direction */
        memcpy(A, A_orig, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
        memcpy(X, X_orig, sizeof(X[0])*(MAX_SIZE));
        /* Random direction for this idir */
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
        for (j = 0; j < n; j++) { for (i = 0; i < n; i++) { if (uplo == CblasUpper && i > j) A_dir[i + j*lda] = 0.0; if (uplo == CblasLower && i < j) A_dir[i + j*lda] = 0.0; } if (diag == CblasUnit) A_dir[j + j*lda] = 0.0; }
        for (i = 0; i < MAX_SIZE; i++) X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
        /* Forward */
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] + h * A_dir[i];
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
        cblas_ztrsv(layout, uplo, transa, diag, n, A, lda, X, incX);
        memcpy(X_plus, X, sizeof(X[0])*(MAX_SIZE));
        /* Backward */
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
        cblas_ztrsv(layout, uplo, transa, diag, n, A, lda, X, incX);
        memcpy(X_minus, X, sizeof(X[0])*(MAX_SIZE));

        vjp_fd = 0.0;
        {
            double temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(X_b_orig[i][idir]) * ((X_plus[i] - X_minus[i]) / (2.0*h)));
            qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
            for (idx = 0; idx < MAX_SIZE; idx++) vjp_fd += temp_products[idx];
        }
        vjp_ad = 0.0;
        {
            double temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = creal(conj(A_dir[i]) * A_b[i][idir]);
            qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(double), compare_abs_d);
            for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }
        {
            double temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(X_dir[i]) * X_b[i][idir]);
            qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
            for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
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

