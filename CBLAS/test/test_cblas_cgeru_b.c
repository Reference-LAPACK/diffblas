/* Test program for cblas_cgeru reverse mode (VJP verification, generic) */
/* Generated automatically by run_tapenade_cblas.py - same methodology as BLAS test_*_reverse.f90 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "cblas.h"
#include "cblas_f77.h"

#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE
extern void set_isize1ofx_(int *val);
extern void set_isize1ofy_(int *val);

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_cgeru(const CBLAS_LAYOUT layout, const CBLAS_INT M, const CBLAS_INT N, const void *alpha, const void *X, const CBLAS_INT incX, const void *Y, const CBLAS_INT incY, void *A, const CBLAS_INT lda);
extern void cblas_cgeru_b(const CBLAS_LAYOUT layout, const CBLAS_INT M, const CBLAS_INT N, const void *alpha, void *alpha_b, const void *X, void *X_b, const CBLAS_INT incX, const void *Y, void *Y_b, const CBLAS_INT incY, void *A, void *A_b, const CBLAS_INT lda);

int main(void) {
    int i, j, idx, n_products;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofx_(&diffblas_isize);
        set_isize1ofy_(&diffblas_isize);
    }
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;
    float vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_INT m = TEST_SIZE;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT incX = 1;
    CBLAS_INT incY = 1;

    float complex alpha[MAX_SIZE], alpha_b[MAX_SIZE], alpha_orig[MAX_SIZE], alpha_dir[MAX_SIZE];
    float complex X[MAX_SIZE], X_b[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    float complex Y[MAX_SIZE], Y_b[MAX_SIZE], Y_orig[MAX_SIZE], Y_dir[MAX_SIZE];
    float complex A[MAX_SIZE*MAX_SIZE], A_b[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    float complex A_plus[MAX_SIZE*MAX_SIZE], A_minus[MAX_SIZE*MAX_SIZE], A_central_diff[MAX_SIZE*MAX_SIZE], A_b_orig[MAX_SIZE*MAX_SIZE];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) { alpha[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); alpha_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { Y[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }

    memcpy(alpha_orig, alpha, sizeof(alpha[0])*(MAX_SIZE));
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));
    memcpy(Y_orig, Y, sizeof(Y[0])*(MAX_SIZE));
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));

    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); A_b_orig[i] = A_b[i]; }
    for (i = 0; i < MAX_SIZE; i++) alpha_b[i] = 0.0f;
    for (i = 0; i < MAX_SIZE; i++) X_b[i] = 0.0f;
    for (i = 0; i < MAX_SIZE; i++) Y_b[i] = 0.0f;

    cblas_cgeru_b(layout, m, n, &alpha, alpha_b, X, X_b, incX, Y, Y_b, incY, A, A_b, lda);

    for (i = 0; i < MAX_SIZE; i++) alpha[i] = alpha_orig[i] + h * alpha_dir[i];
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
    for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] + h * Y_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] + h * A_dir[i];
    cblas_cgeru(layout, m, n, &alpha, X, incX, Y, incY, A, lda);
    memcpy(A_plus, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) alpha[i] = alpha_orig[i] - h * alpha_dir[i];
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
    for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] - h * Y_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
    cblas_cgeru(layout, m, n, &alpha, X, incX, Y, incY, A, lda);
    memcpy(A_minus, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));

    vjp_fd = 0.0f;
    {
        float temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = creal(conj(A_b_orig[i]) * ((A_plus[i] - A_minus[i]) / (2.0*h)));
        qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_fd += temp_products[idx];
    }

    vjp_ad = 0.0f;
    {
        float temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(alpha_dir[i]) * alpha_b[i]);
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        float temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(X_dir[i]) * X_b[i]);
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        float temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(Y_dir[i]) * Y_b[i]);
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        float temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = creal(conj(A_dir[i]) * A_b[i]);
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

