/* Test program for cblas_ssyr2k vector reverse mode (VJP verification, generic, loop over directions) */
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
extern void set_isize2ofa_(int *val);
extern void set_isize2ofb_(int *val);

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_ssyr2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE Trans, const CBLAS_INT N, const CBLAS_INT K, const float alpha, const float *A, const CBLAS_INT lda, const float *B, const CBLAS_INT ldb, const float beta, float *C, const CBLAS_INT ldc);
/* cblas_*_bv from cblas_bv.h */

int main(void) {
    int i, j, idx, idir, nbdirs = NBDirsMax, n_products;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize2ofa_(&diffblas_isize);
        set_isize2ofb_(&diffblas_isize);
    }
    int has_large_errors = 0;
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;
    float max_error = 0.0f;
    float vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT k = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT ldb = MAX_SIZE;
    CBLAS_INT ldc = MAX_SIZE;

    float alpha, alpha_b[NBDirsMax], alpha_orig, alpha_dir, alpha_b_orig[NBDirsMax];
    float A[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    float A_b[MAX_SIZE*MAX_SIZE][NBDirsMax], A_b_orig[MAX_SIZE*MAX_SIZE][NBDirsMax];
    float B[MAX_SIZE*MAX_SIZE], B_orig[MAX_SIZE*MAX_SIZE], B_dir[MAX_SIZE*MAX_SIZE];
    float B_b[MAX_SIZE*MAX_SIZE][NBDirsMax], B_b_orig[MAX_SIZE*MAX_SIZE][NBDirsMax];
    float beta, beta_b[NBDirsMax], beta_orig, beta_dir, beta_b_orig[NBDirsMax];
    float C[MAX_SIZE*MAX_SIZE], C_orig[MAX_SIZE*MAX_SIZE], C_dir[MAX_SIZE*MAX_SIZE];
    float C_b[MAX_SIZE*MAX_SIZE][NBDirsMax], C_b_orig[MAX_SIZE*MAX_SIZE][NBDirsMax];
    float C_plus[MAX_SIZE*MAX_SIZE], C_minus[MAX_SIZE*MAX_SIZE];

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
    for (j = 0; j < NBDirsMax; j++) alpha_b[j] = 0.0f;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) A_b[i][j] = 0.0f;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) B_b[i][j] = 0.0f;
    for (j = 0; j < NBDirsMax; j++) beta_b[j] = 0.0f;

    cblas_ssyr2k_bv(layout, uplo, trans, n, k, alpha, &alpha_b, A, &A_b[0][0], lda, B, &B_b[0][0], ldb, beta, &beta_b, C, &C_b[0][0], ldc, nbdirs);

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
        cblas_ssyr2k(layout, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        memcpy(C_plus, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));
        /* Backward */
        alpha = alpha_orig - h * alpha_dir;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) B[i] = B_orig[i] - h * B_dir[i];
        beta = beta_orig - h * beta_dir;
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) C[i] = C_orig[i] - h * C_dir[i];
        cblas_ssyr2k(layout, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        memcpy(C_minus, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));

        vjp_fd = 0.0f;
        {
            float temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = C_b_orig[i][idir] * ((C_plus[i] - C_minus[i]) / (2.0*h));
            qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_fd += temp_products[idx];
        }
        vjp_ad = 0.0f;
        vjp_ad += alpha_dir * alpha_b[idir];
        {
            float temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = A_dir[i] * A_b[i][idir];
            qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }
        {
            float temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = B_dir[i] * B_b[i][idir];
            qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }
        vjp_ad += beta_dir * beta_b[idir];
        {
            float temp_products[MAX_SIZE*MAX_SIZE];
            for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = C_dir[i] * C_b[i][idir];
            qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
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

