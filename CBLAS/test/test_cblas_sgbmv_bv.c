/* Test program for cblas_sgbmv vector reverse mode (VJP verification, generic, loop over directions) */
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
extern void set_isize1ofx_(int *val);
extern void set_isize2ofa_(int *val);

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_sgbmv(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA, const CBLAS_INT M, const CBLAS_INT N, const CBLAS_INT KL, const CBLAS_INT KU, const float alpha, const float *A, const CBLAS_INT lda, const float *X, const CBLAS_INT incX, const float beta, float *Y, const CBLAS_INT incY);
/* cblas_*_bv from cblas_bv.h */

int main(void) {
    int i, j, idx, idir, nbdirs = NBDirsMax, n_products;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofx_(&diffblas_isize);
        set_isize2ofa_(&diffblas_isize);
    }
    int has_large_errors = 0;
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;
    float max_error = 0.0f;
    float vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE transa = CblasNoTrans;
    CBLAS_INT m = TEST_SIZE;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT incX = 1;
    CBLAS_INT incY = 1;

    CBLAS_INT KL = 1;  /* band width: LDA >= KL+KU+1 (match _d/dv/b) */
    CBLAS_INT KU = 1;  /* band width: LDA >= KL+KU+1 (match _d/dv/b) */
    float alpha, alpha_b[NBDirsMax], alpha_orig, alpha_dir, alpha_b_orig[NBDirsMax];
    float A[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    float A_b[MAX_SIZE*MAX_SIZE][NBDirsMax], A_b_orig[MAX_SIZE*MAX_SIZE][NBDirsMax];
    float X[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    float X_b[MAX_SIZE][NBDirsMax], X_b_orig[MAX_SIZE][NBDirsMax];
    float beta, beta_b[NBDirsMax], beta_orig, beta_dir, beta_b_orig[NBDirsMax];
    float Y[MAX_SIZE], Y_orig[MAX_SIZE], Y_dir[MAX_SIZE];
    float Y_b[MAX_SIZE][NBDirsMax], Y_b_orig[MAX_SIZE][NBDirsMax];
    float Y_plus[MAX_SIZE], Y_minus[MAX_SIZE];

    srand(42);
    alpha = (rand()/(double)RAND_MAX)*2.0 - 1.0; alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    /* A: general band storage (KL+KU+1) x N (match BLAS/test) */
    memset(A, 0, sizeof(A));
    for (j = 0; j < MAX_SIZE; j++) {
        int band_rows = KL + KU + 1;
        for (i = 0; i < band_rows; i++) {
            A[i + j * MAX_SIZE] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
        }
    }
    /* A: general band storage (KL+KU+1) x N (match BLAS/test) */
    memset(A_dir, 0, sizeof(A_dir));
    for (j = 0; j < MAX_SIZE; j++) {
        int band_rows = KL + KU + 1;
        for (i = 0; i < band_rows; i++) {
            A_dir[i + j * MAX_SIZE] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
        }
    }
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    beta = (rand()/(double)RAND_MAX)*2.0 - 1.0; beta_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) { Y[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }

    alpha_orig = alpha;
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));
    beta_orig = beta;
    memcpy(Y_orig, Y, sizeof(Y[0])*(MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) { Y_b[i][j] = (rand()/(double)RAND_MAX)*2.0 - 1.0; Y_b_orig[i][j] = Y_b[i][j]; }
    for (j = 0; j < NBDirsMax; j++) alpha_b[j] = 0.0f;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) A_b[i][j] = 0.0f;
    for (i = 0; i < MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) X_b[i][j] = 0.0f;
    for (j = 0; j < NBDirsMax; j++) beta_b[j] = 0.0f;

    cblas_sgbmv_bv(layout, transa, m, n, KL, KU, alpha, &alpha_b, A, &A_b[0][0], lda, X, X_b, incX, beta, &beta_b, Y, Y_b, incY, nbdirs);

    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals for this direction */
        alpha = alpha_orig;
        memcpy(A, A_orig, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
        memcpy(X, X_orig, sizeof(X[0])*(MAX_SIZE));
        beta = beta_orig;
        memcpy(Y, Y_orig, sizeof(Y[0])*(MAX_SIZE));
        /* Random direction for this idir */
        alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        /* A: general band storage (KL+KU+1) x N (match BLAS/test) */
        memset(A_dir, 0, sizeof(A_dir));
        for (j = 0; j < MAX_SIZE; j++) {
            int band_rows = KL + KU + 1;
            for (i = 0; i < band_rows; i++) {
                A_dir[i + j * MAX_SIZE] = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
            }
        }
        for (i = 0; i < MAX_SIZE; i++) X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        beta_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        for (i = 0; i < MAX_SIZE; i++) Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        /* Forward */
        alpha = alpha_orig + h * alpha_dir;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] + h * A_dir[i];
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
        beta = beta_orig + h * beta_dir;
        for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] + h * Y_dir[i];
        cblas_sgbmv(layout, transa, m, n, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
        memcpy(Y_plus, Y, sizeof(Y[0])*(MAX_SIZE));
        /* Backward */
        alpha = alpha_orig - h * alpha_dir;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
        beta = beta_orig - h * beta_dir;
        for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] - h * Y_dir[i];
        cblas_sgbmv(layout, transa, m, n, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
        memcpy(Y_minus, Y, sizeof(Y[0])*(MAX_SIZE));

        vjp_fd = 0.0f;
        {
            float temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = Y_b_orig[i][idir] * ((Y_plus[i] - Y_minus[i]) / (2.0*h));
            qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < MAX_SIZE; idx++) vjp_fd += temp_products[idx];
        }
        vjp_ad = 0.0f;
        vjp_ad += alpha_dir * alpha_b[idir];
        {
            float temp_products[MAX_SIZE*MAX_SIZE];
            int n_band = 0;
            int band_rows = KL + KU + 1;
            for (j = 0; j < n; j++) for (i = 0; i < band_rows; i++) {
                temp_products[n_band++] = A_dir[i+j*lda] * A_b[i+j*lda][idir];
            }
            qsort(temp_products, (size_t)n_band, sizeof(float), compare_abs_f);
            for (idx = 0; idx < n_band; idx++) vjp_ad += temp_products[idx];
        }
        {
            float temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = X_dir[i] * X_b[i][idir];
            qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }
        vjp_ad += beta_dir * beta_b[idir];
        {
            float temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = Y_dir[i] * Y_b[i][idir];
            qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }

        {
            float abs_err = fabsf(vjp_fd - vjp_ad);
            float abs_reference = fabsf(vjp_ad);
            float error_bound = atol + rtol * (abs_reference > 1e-10f ? abs_reference : 1e-10f);
            if (abs_err > error_bound) has_large_errors = 1;
            { float r = abs_err / error_bound; if (r > max_error) max_error = r; }
        }
    }

    printf("Maximum error ratio: %.6e\n", (double)max_error);
    if (has_large_errors) { printf("FAIL: Large errors in derivatives\n"); return 1; }
    if (max_error < 0.5f) { printf("PASS: Derivatives accurate to machine precision\n"); return 0; }
    if (max_error < 1.0f) { printf("PASS: Derivatives reasonably accurate\n"); return 0; }
    printf("WARNING: Derivatives may have significant errors\n"); return 0;
}

