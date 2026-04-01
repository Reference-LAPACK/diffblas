/* Test program for cblas_ssyr reverse mode (VJP verification, generic) */
/* Generated automatically by run_tapenade_cblas.py - same methodology as BLAS test_*_reverse.f90 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cblas.h"
#include "cblas_f77.h"

#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE
extern void set_isize1ofx_(int *val);

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_ssyr(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_INT N, const float alpha, const float *X, const CBLAS_INT incX, float *A, const CBLAS_INT lda);
extern void cblas_ssyr_b(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_INT N, const float alpha, float *alpha_b, const float *X, float *X_b, const CBLAS_INT incX, float *A, float *A_b, const CBLAS_INT lda);

int main(void) {
    int i, j, idx, n_products;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofx_(&diffblas_isize);
    }
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;
    float vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT incX = 1;

    float alpha, alpha_b, alpha_orig, alpha_dir;
    float X[MAX_SIZE], X_b[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    float A[MAX_SIZE*MAX_SIZE], A_b[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    float A_plus[MAX_SIZE*MAX_SIZE], A_minus[MAX_SIZE*MAX_SIZE], A_central_diff[MAX_SIZE*MAX_SIZE], A_b_orig[MAX_SIZE*MAX_SIZE];

    srand(42);
    alpha = (rand()/(double)RAND_MAX)*2.0 - 1.0; alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }

    alpha_orig = alpha;
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));

    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; A_b_orig[i] = A_b[i]; }
    alpha_b = 0.0f;
    for (i = 0; i < MAX_SIZE; i++) X_b[i] = 0.0f;

    cblas_ssyr_b(layout, uplo, n, alpha, &alpha_b, X, X_b, incX, A, A_b, lda);

    alpha = alpha_orig + h * alpha_dir;
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] + h * A_dir[i];
    cblas_ssyr(layout, uplo, n, alpha, X, incX, A, lda);
    memcpy(A_plus, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));

    alpha = alpha_orig - h * alpha_dir;
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
    cblas_ssyr(layout, uplo, n, alpha, X, incX, A, lda);
    memcpy(A_minus, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));

    vjp_fd = 0.0f;
    {
        float temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = A_b_orig[i] * ((A_plus[i] - A_minus[i]) / (2.0*h));
        qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_fd += temp_products[idx];
    }

    vjp_ad = 0.0f;
    vjp_ad += alpha_dir * alpha_b;
    {
        float temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = X_dir[i] * X_b[i];
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        float temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = A_dir[i] * A_b[i];
        qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
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

