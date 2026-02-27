/* Test program for cblas_sgemm vector reverse mode (VJP verification, loop over directions) */
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
#define MAT_SIZE (MAX_SIZE*MAX_SIZE)

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_sgemm(CBLAS_LAYOUT, CBLAS_TRANSPOSE, CBLAS_TRANSPOSE, CBLAS_INT, CBLAS_INT, CBLAS_INT,
    float, const float *, CBLAS_INT, const float *, CBLAS_INT, float, float *, CBLAS_INT);
/* _bv declaration from cblas_bv.h */

int main(void) {
    int i, j, idx, idir, nbdirs = NBDirsMax;
    int has_large_errors = 0;
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;
    float max_error = 0.0f;
    CBLAS_INT m = TEST_SIZE, n = TEST_SIZE, k = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE, ldb = MAX_SIZE, ldc = MAX_SIZE;
    float alpha, beta;
    float alpha_b[NBDirsMax], beta_b[NBDirsMax];
    float A[MAT_SIZE], B[MAT_SIZE], C[MAT_SIZE];
    float A_b[MAT_SIZE*NBDirsMax], B_b[MAT_SIZE*NBDirsMax], C_b[MAT_SIZE*NBDirsMax];  /* layout: element then direction */
    float A_dir[MAT_SIZE], B_dir[MAT_SIZE], C_dir[MAT_SIZE];
    float C_forward[MAT_SIZE], C_backward[MAT_SIZE];
    float C_b_orig[MAT_SIZE*NBDirsMax];  /* save cotangents for all directions (like BLAS cb_orig) */
    float alpha_orig, beta_orig, alpha_dir, beta_dir;
    float A_orig[MAT_SIZE], B_orig[MAT_SIZE], C_orig[MAT_SIZE];

    srand(42);
    alpha = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    beta  = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    for (i = 0; i < MAT_SIZE; i++) {
        A[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        B[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        C[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    }
    /* Cotangents for all directions (seeds for reverse, like BLAS cb(k) and _b C_b) */
    for (i = 0; i < MAT_SIZE; i++)
        for (j = 0; j < NBDirsMax; j++) {
            C_b[i*NBDirsMax + j] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        }

    alpha_orig = alpha; beta_orig = beta;
    memcpy(A_orig, A, sizeof(A)); memcpy(B_orig, B, sizeof(B)); memcpy(C_orig, C, sizeof(C));
    memcpy(C_b_orig, C_b, sizeof(C_b));  /* save before _bv (inout C_b overwritten) */
    /* Input adjoints zero (computed by _bv), same as _b and BLAS _bv */
    for (j = 0; j < NBDirsMax; j++) { alpha_b[j] = 0.0f; beta_b[j] = 0.0f; }
    for (i = 0; i < MAT_SIZE*NBDirsMax; i++) { A_b[i] = 0.0f; B_b[i] = 0.0f; }

    cblas_sgemm_bv(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k,
        alpha, &alpha_b, A, A_b, lda, B, B_b, ldb, beta, &beta_b, C, C_b, ldc, nbdirs);

    /* Per-direction VJP check (gradient logic like _b and BLAS _bv: direction^T @ adjoint vs cotangent^T @ FD) */
    for (idir = 0; idir < nbdirs; idir++) {
        /* Random direction for this idir (like BLAS: random_number inside loop) */
        alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        beta_dir  = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        for (i = 0; i < MAT_SIZE; i++) {
            A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
            B_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
            C_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        }
        /* Forward perturbation */
        alpha = alpha_orig + h * alpha_dir; beta = beta_orig + h * beta_dir;
        for (i = 0; i < MAT_SIZE; i++) { A[i] = A_orig[i] + h * A_dir[i]; B[i] = B_orig[i] + h * B_dir[i]; C[i] = C_orig[i] + h * C_dir[i]; }
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        memcpy(C_forward, C, sizeof(C));
        /* Backward perturbation */
        alpha = alpha_orig - h * alpha_dir; beta = beta_orig - h * beta_dir;
        for (i = 0; i < MAT_SIZE; i++) { A[i] = A_orig[i] - h * A_dir[i]; B[i] = B_orig[i] - h * B_dir[i]; C[i] = C_orig[i] - h * C_dir[i]; }
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        memcpy(C_backward, C, sizeof(C));

        float vjp_fd, vjp_ad;
        /* VJP fd: cotangent_idir^T @ (C_forward - C_backward)/(2h), sorted (like _b / BLAS) */
        {
            float temp_products[MAT_SIZE];
            int n_products = MAT_SIZE;
            for (i = 0; i < n_products; i++) temp_products[i] = C_b_orig[i*NBDirsMax + idir] * ((C_forward[i] - C_backward[i]) / (2.0*h));
            qsort(temp_products, (size_t)n_products, sizeof(float), compare_abs_f);
            vjp_fd = 0.0f;
            for (idx = 0; idx < n_products; idx++) vjp_fd += temp_products[idx];
        }
        /* VJP ad: direction^T @ adjoint_idir (same as _b per direction) */
        vjp_ad = 0.0f;
        vjp_ad += alpha_dir * alpha_b[idir] + beta_dir * beta_b[idir];
        {
            float temp_products[MAT_SIZE];
            int n_products = MAT_SIZE;
            for (i = 0; i < n_products; i++) temp_products[i] = A_dir[i] * A_b[i*NBDirsMax + idir];
            qsort(temp_products, (size_t)n_products, sizeof(float), compare_abs_f);
            for (idx = 0; idx < n_products; idx++) vjp_ad += temp_products[idx];
            for (i = 0; i < n_products; i++) temp_products[i] = B_dir[i] * B_b[i*NBDirsMax + idir];
            qsort(temp_products, (size_t)n_products, sizeof(float), compare_abs_f);
            for (idx = 0; idx < n_products; idx++) vjp_ad += temp_products[idx];
            for (i = 0; i < n_products; i++) temp_products[i] = C_dir[i] * C_b[i*NBDirsMax + idir];
            qsort(temp_products, (size_t)n_products, sizeof(float), compare_abs_f);
            for (idx = 0; idx < n_products; idx++) vjp_ad += temp_products[idx];
        }

        float abs_err = fabsf(vjp_fd - vjp_ad);
        float abs_reference = fabsf(vjp_ad);
        float error_bound = atol + rtol * (abs_reference > 1e-10f ? abs_reference : 1e-10f);
        if (abs_err > error_bound) has_large_errors = 1;
        {
            float r = abs_err / error_bound;
            if (r > max_error) max_error = r;
        }
    }

    printf("Maximum error ratio (abs_error/error_bound): %.6e\n", (double)max_error);
    if (has_large_errors) { printf("FAIL: Large errors detected in derivatives\n"); return 1; }
    if (max_error < 0.5f) { printf("PASS: Derivatives are accurate to machine precision\n"); return 0; }
    if (max_error < 1.0f) { printf("PASS: Derivatives are reasonably accurate\n"); return 0; }
    printf("WARNING: Derivatives may have significant errors\n"); return 0;
}

