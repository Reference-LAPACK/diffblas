/* Test program for cblas_dasum reverse mode (VJP verification, generic) */
/* Generated automatically by run_tapenade_cblas.py - same methodology as BLAS test_*_reverse.f90 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cblas.h"
#include "cblas_f77.h"

#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE

static int compare_abs_d(const void *a, const void *b) { double x = fabs(*(const double*)a), y = fabs(*(const double*)b); return (x > y) - (x < y); }

extern double cblas_dasum(const CBLAS_INT N, const double *X, const CBLAS_INT incX);
extern void cblas_dasum_b(const CBLAS_INT N, const double *X, double *X_b, const CBLAS_INT incX);

int main(void) {
    int i, j, idx, n_products;
    double h = 1.0e-7;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double vjp_fd, vjp_ad;

    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT incX = 1;

    double X[MAX_SIZE], X_b[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }

    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) X_b[i] = 0.0;

    cblas_dasum_b(n, X, X_b, incX);

    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
    cblas_dasum(n, X, incX);

    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
    cblas_dasum(n, X, incX);

    vjp_fd = 0.0;

    vjp_ad = 0.0;
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = X_dir[i] * X_b[i];
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

