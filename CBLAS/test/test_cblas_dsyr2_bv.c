/* Test program for cblas_dsyr2 vector reverse mode (VJP verification, generic, loop over directions) */
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

extern void cblas_dsyr2(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_INT N, const double alpha, const double *X, const CBLAS_INT incX, const double *Y, const CBLAS_INT incY, double *A, const CBLAS_INT lda);
/* cblas_*_bv from cblas_bv.h */

int main(void) {
    int i, j, idx, idir, nbdirs = NBDirsMax, n_products;
    int has_large_errors = 0;
    double h = 1.0e-7;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double max_error = 0.0;
    double vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT incX = 1;
    CBLAS_INT incY = 1;

    double alpha, alpha_b[NBDirsMax], alpha_orig, alpha_dir, alpha_b_orig[NBDirsMax];
    double X[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    double X_b[MAX_SIZE][NBDirsMax], X_b_orig[MAX_SIZE][NBDirsMax];
    double Y[MAX_SIZE], Y_orig[MAX_SIZE], Y_dir[MAX_SIZE];
    double Y_b[MAX_SIZE][NBDirsMax], Y_b_orig[MAX_SIZE][NBDirsMax];
    double A[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    double A_b[MAX_SIZE*MAX_SIZE][NBDirsMax], A_b_orig[MAX_SIZE*MAX_SIZE][NBDirsMax];
    double A_plus[MAX_SIZE*MAX_SIZE], A_minus[MAX_SIZE*MAX_SIZE];

    srand(42);
    alpha = (rand()/(double)RAND_MAX)*2.0 - 1.0; alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE; i++) { Y[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }

    alpha_orig = alpha;
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));
    memcpy(Y_orig, Y, sizeof(Y[0])*(MAX_SIZE));
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));

    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) { A_b[i][j] = (rand()/(double)RAND_MAX)*2.0 - 1.0; A_b_orig[i][j] = A_b[i][j]; }
    for (j = 0; j < NBDirsMax; j++) alpha_b[j] = 0.0;
    for (i = 0; i < MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) X_b[i][j] = 0.0;
    for (i = 0; i < MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) Y_b[i][j] = 0.0;

    cblas_dsyr2_bv(layout, uplo, n, alpha, &alpha_b, X, X_b, incX, Y, Y_b, incY, A, A_b, lda, nbdirs);

    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals for this direction */
        alpha = alpha_orig;
        memcpy(X, X_orig, sizeof(X[0])*(MAX_SIZE));
        memcpy(Y, Y_orig, sizeof(Y[0])*(MAX_SIZE));
        memcpy(A, A_orig, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
        /* Random direction for this idir */
        alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        for (i = 0; i < MAX_SIZE; i++) X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        for (i = 0; i < MAX_SIZE; i++) Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        /* Forward */
        alpha = alpha_orig + h * alpha_dir;
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
        for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] + h * Y_dir[i];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] + h * A_dir[i];
        cblas_dsyr2(layout, uplo, n, alpha, X, incX, Y, incY, A, lda);
        memcpy(A_plus, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
        /* Backward */
        alpha = alpha_orig - h * alpha_dir;
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
        for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] - h * Y_dir[i];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
        cblas_dsyr2(layout, uplo, n, alpha, X, incX, Y, incY, A, lda);
        memcpy(A_minus, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));

        vjp_fd = 0.0;
        {
            double temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = A_b_orig[i][idir] * ((A_plus[i] - A_minus[i]) / (2.0*h));
            qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(double), compare_abs_d);
            for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_fd += temp_products[idx];
        }
        vjp_ad = 0.0;
        vjp_ad += alpha_dir * alpha_b[idir];
        {
            double temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = X_dir[i] * X_b[i][idir];
            qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
            for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }
        {
            double temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = Y_dir[i] * Y_b[i][idir];
            qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
            for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }
        {
            double temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = A_dir[i] * A_b[i][idir];
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

