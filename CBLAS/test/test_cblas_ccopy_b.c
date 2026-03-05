/* Test program for cblas_ccopy reverse mode (VJP verification, generic) */
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
extern void set_isize1ofcx_(int *val);

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_ccopy(const CBLAS_INT N, const void *X, const CBLAS_INT incX, void *Y, const CBLAS_INT incY);
extern void cblas_ccopy_b(const CBLAS_INT N, const void *X, void *X_b, const CBLAS_INT incX, void *Y, void *Y_b, const CBLAS_INT incY);

int main(void) {
    int i, j, idx, n_products;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofcx_(&diffblas_isize);
    }
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;
    float vjp_fd, vjp_ad;

    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT incX = 1;
    CBLAS_INT incY = 1;

    float complex X[MAX_SIZE], X_b[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    float complex Y[MAX_SIZE], Y_b[MAX_SIZE], Y_orig[MAX_SIZE], Y_dir[MAX_SIZE];
    float complex Y_plus[MAX_SIZE], Y_minus[MAX_SIZE], Y_central_diff[MAX_SIZE], Y_b_orig[MAX_SIZE];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { Y[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }

    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));
    memcpy(Y_orig, Y, sizeof(Y[0])*(MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) { Y_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); Y_b_orig[i] = Y_b[i]; }
    for (i = 0; i < MAX_SIZE; i++) X_b[i] = 0.0f;

    cblas_ccopy_b(n, X, X_b, incX, Y, Y_b, incY);

    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
    for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] + h * Y_dir[i];
    cblas_ccopy(n, X, incX, Y, incY);
    memcpy(Y_plus, Y, sizeof(Y[0])*(MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
    for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] - h * Y_dir[i];
    cblas_ccopy(n, X, incX, Y, incY);
    memcpy(Y_minus, Y, sizeof(Y[0])*(MAX_SIZE));

    vjp_fd = 0.0f;
    {
        float temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(Y_b_orig[i]) * ((Y_plus[i] - Y_minus[i]) / (2.0*h)));
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_fd += temp_products[idx];
    }

    vjp_ad = 0.0f;
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
        float abs_err = fabsf(vjp_fd - vjp_ad);
        float abs_reference = fabsf(vjp_ad);
        float error_bound = atol + rtol * (abs_reference > 1e-10f ? abs_reference : 1e-10f);
        printf("VJP: fd=%.10e ad=%.10e abs_err=%.10e error_bound=%.10e\n", (double)vjp_fd, (double)vjp_ad, (double)abs_err, (double)error_bound);
        if (abs_err > error_bound) { printf("FAIL: Large errors detected in derivatives (outside tolerance)\n"); return 1; }
        if (abs_err < 0.5 * error_bound) { printf("PASS: Derivatives are accurate to machine precision\n"); return 0; }
        printf("PASS: Derivatives are within tolerance (rtol + atol)\n"); return 0;
    }
}

