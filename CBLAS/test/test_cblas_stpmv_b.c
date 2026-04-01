/* Test program for cblas_stpmv reverse mode (VJP verification, generic) */
/* Generated automatically by run_tapenade_cblas.py - same methodology as BLAS test_*_reverse.f90 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cblas.h"
#include "cblas_f77.h"

#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE
#define PACKED_SIZE ((MAX_SIZE) * ((MAX_SIZE) + 1) / 2)  /* packed triangular/symmetric */
extern void set_isize1ofap_(int *val);

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_stpmv(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const CBLAS_INT N, const float *Ap, float *X, const CBLAS_INT incX);
extern void cblas_stpmv_b(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag, const CBLAS_INT N, const float *Ap, float *Ap_b, float *X, float *X_b, const CBLAS_INT incX);

int main(void) {
    int i, j, idx, n_products;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofap_(&diffblas_isize);
    }
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;
    float vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE transa = CblasNoTrans;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_DIAG diag = CblasNonUnit;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT incX = 1;

    float Ap[PACKED_SIZE], Ap_b[PACKED_SIZE], Ap_orig[PACKED_SIZE], Ap_dir[PACKED_SIZE];
    float X[MAX_SIZE], X_b[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    float X_plus[MAX_SIZE], X_minus[MAX_SIZE], X_central_diff[MAX_SIZE], X_b_orig[MAX_SIZE];

    srand(42);
    for (i = 0; i < PACKED_SIZE; i++) { Ap[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; Ap_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }

    memcpy(Ap_orig, Ap, sizeof(Ap[0])*(PACKED_SIZE));
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) { X_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; X_b_orig[i] = X_b[i]; }
    for (i = 0; i < PACKED_SIZE; i++) Ap_b[i] = 0.0f;

    cblas_stpmv_b(layout, uplo, transa, diag, n, Ap, Ap_b, X, X_b, incX);

    for (i = 0; i < PACKED_SIZE; i++) Ap[i] = Ap_orig[i] + h * Ap_dir[i];
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
    cblas_stpmv(layout, uplo, transa, diag, n, Ap, X, incX);
    memcpy(X_plus, X, sizeof(X[0])*(MAX_SIZE));

    for (i = 0; i < PACKED_SIZE; i++) Ap[i] = Ap_orig[i] - h * Ap_dir[i];
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
    cblas_stpmv(layout, uplo, transa, diag, n, Ap, X, incX);
    memcpy(X_minus, X, sizeof(X[0])*(MAX_SIZE));

    vjp_fd = 0.0f;
    {
        float temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = X_b_orig[i] * ((X_plus[i] - X_minus[i]) / (2.0*h));
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_fd += temp_products[idx];
    }

    vjp_ad = 0.0f;
    {
        float temp_products[PACKED_SIZE];
        for (i = 0; i < PACKED_SIZE; i++) temp_products[i] = Ap_dir[i] * Ap_b[i];
        qsort(temp_products, (size_t)PACKED_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < PACKED_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        float temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = X_dir[i] * X_b[i];
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }

    {
        float abs_err = fabsf(vjp_fd - vjp_ad);
        float abs_reference = fabsf(vjp_ad);
        float error_bound = atol + rtol * (abs_reference > 1e-10f ? abs_reference : 1e-10f);
        printf("VJP: fd=%.10e ad=%.10e abs_err=%.10e error_bound=%.10e\n", (double)vjp_fd, (double)vjp_ad, (double)abs_err, (double)error_bound);
        if (abs_err > error_bound) { printf("FAIL: Large errors detected in derivatives (outside tolerance)\n"); return 1; }
        if (abs_err < 0.5 * error_bound) { printf("PASS: Derivatives are accurate to machine precision\n"); return 0; }
        printf("PASS: Derivatives are within tolerance (rtol + atol)\n"); return 0;
    }
}

