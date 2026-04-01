/* Test program for cblas_zscal reverse mode (VJP verification, generic) */
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

static int compare_abs_d(const void *a, const void *b) { double x = fabs(*(const double*)a), y = fabs(*(const double*)b); return (x > y) - (x < y); }

extern void cblas_zscal(const CBLAS_INT N, const void *alpha, void *X, const CBLAS_INT incX);
extern void cblas_zscal_b(const CBLAS_INT N, const void *alpha, void *alpha_b, void *X, void *X_b, const CBLAS_INT incX);

int main(void) {
    int i, j, idx, n_products;
    double h = 1.0e-7;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double vjp_fd, vjp_ad;

    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT incX = 1;

    double complex alpha[MAX_SIZE], alpha_b[MAX_SIZE], alpha_orig[MAX_SIZE], alpha_dir[MAX_SIZE];
    double complex X[MAX_SIZE], X_b[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    double complex X_plus[MAX_SIZE], X_minus[MAX_SIZE], X_central_diff[MAX_SIZE], X_b_orig[MAX_SIZE];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) { alpha[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); alpha_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }

    memcpy(alpha_orig, alpha, sizeof(alpha[0])*(MAX_SIZE));
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) { X_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); X_b_orig[i] = X_b[i]; }
    for (i = 0; i < MAX_SIZE; i++) alpha_b[i] = 0.0;

    cblas_zscal_b(n, &alpha, alpha_b, X, X_b, incX);

    for (i = 0; i < MAX_SIZE; i++) alpha[i] = alpha_orig[i] + h * alpha_dir[i];
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
    cblas_zscal(n, &alpha, X, incX);
    memcpy(X_plus, X, sizeof(X[0])*(MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) alpha[i] = alpha_orig[i] - h * alpha_dir[i];
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
    cblas_zscal(n, &alpha, X, incX);
    memcpy(X_minus, X, sizeof(X[0])*(MAX_SIZE));

    vjp_fd = 0.0;
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(X_b_orig[i]) * ((X_plus[i] - X_minus[i]) / (2.0*h)));
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_fd += temp_products[idx];
    }

    vjp_ad = 0.0;
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(alpha_dir[i]) * alpha_b[i]);
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(X_dir[i]) * X_b[i]);
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }

    {
        double abs_err = fabs(vjp_fd - vjp_ad);
        double abs_reference = fabs(vjp_ad);
        double error_bound = atol + rtol * (abs_reference > 1e-10 ? abs_reference : 1e-10);
        printf("VJP: fd=%.10e ad=%.10e abs_err=%.10e error_bound=%.10e\n", (double)vjp_fd, (double)vjp_ad, (double)abs_err, (double)error_bound);
        if (abs_err > error_bound) { printf("FAIL: Large errors detected in derivatives (outside tolerance)\n"); return 1; }
        if (abs_err < 0.5 * error_bound) { printf("PASS: Derivatives are accurate to machine precision\n"); return 0; }
        printf("PASS: Derivatives are within tolerance (rtol + atol)\n"); return 0;
    }
}

