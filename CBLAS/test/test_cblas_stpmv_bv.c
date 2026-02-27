/* Test program for cblas_stpmv vector reverse mode (VJP verification, generic, loop over directions) */
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
#define PACKED_SIZE ((MAX_SIZE) * ((MAX_SIZE) + 1) / 2)

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_stpmv(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const CBLAS_INT N, const float *Ap, float *X, const CBLAS_INT incX);
/* cblas_*_bv from cblas_bv.h */

int main(void) {
    int i, j, idx, idir, nbdirs = NBDirsMax, n_products;
    int has_large_errors = 0;
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;
    float max_error = 0.0f;
    float vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE transa = CblasNoTrans;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_DIAG diag = CblasNonUnit;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT incX = 1;

    float Ap[PACKED_SIZE], Ap_orig[PACKED_SIZE], Ap_dir[PACKED_SIZE];
    float Ap_b[PACKED_SIZE][NBDirsMax], Ap_b_orig[PACKED_SIZE][NBDirsMax];
    float X[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    float X_b[MAX_SIZE][NBDirsMax], X_b_orig[MAX_SIZE][NBDirsMax];
    float X_plus[MAX_SIZE], X_minus[MAX_SIZE];

    srand(42);
    for (i = 0; i < PACKED_SIZE; i++) { Ap[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; Ap_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }

    memcpy(Ap_orig, Ap, sizeof(Ap[0])*(PACKED_SIZE));
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) { X_b[i][j] = (rand()/(double)RAND_MAX)*2.0 - 1.0; X_b_orig[i][j] = X_b[i][j]; }
    for (i = 0; i < PACKED_SIZE; i++) for (j = 0; j < NBDirsMax; j++) Ap_b[i][j] = 0.0f;

    cblas_stpmv_bv(layout, uplo, transa, diag, n, Ap, Ap_b, X, X_b, incX, nbdirs);

    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals for this direction */
        memcpy(Ap, Ap_orig, sizeof(Ap[0])*(PACKED_SIZE));
        memcpy(X, X_orig, sizeof(X[0])*(MAX_SIZE));
        /* Random direction for this idir */
        for (i = 0; i < PACKED_SIZE; i++) Ap_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        for (i = 0; i < MAX_SIZE; i++) X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        /* Forward */
        for (i = 0; i < PACKED_SIZE; i++) Ap[i] = Ap_orig[i] + h * Ap_dir[i];
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
        cblas_stpmv(layout, uplo, transa, diag, n, Ap, X, incX);
        memcpy(X_plus, X, sizeof(X[0])*(MAX_SIZE));
        /* Backward */
        for (i = 0; i < PACKED_SIZE; i++) Ap[i] = Ap_orig[i] - h * Ap_dir[i];
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
        cblas_stpmv(layout, uplo, transa, diag, n, Ap, X, incX);
        memcpy(X_minus, X, sizeof(X[0])*(MAX_SIZE));

        vjp_fd = 0.0f;
        {
            float temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = X_b_orig[i][idir] * ((X_plus[i] - X_minus[i]) / (2.0*h));
            qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < MAX_SIZE; idx++) vjp_fd += temp_products[idx];
        }
        vjp_ad = 0.0f;
        {
            float temp_products[PACKED_SIZE];
            for (i = 0; i < PACKED_SIZE; i++) temp_products[i] = Ap_dir[i] * Ap_b[i][idir];
            qsort(temp_products, (size_t)PACKED_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < PACKED_SIZE; idx++) vjp_ad += temp_products[idx];
        }
        {
            float temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = X_dir[i] * X_b[i][idir];
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

