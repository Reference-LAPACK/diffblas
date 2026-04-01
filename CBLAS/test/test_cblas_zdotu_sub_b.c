/* Test program for cblas_zdotu_sub reverse mode (VJP verification, generic) */
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
extern void set_isize1ofzx_(int *val);
extern void set_isize1ofzy_(int *val);

static int compare_abs_d(const void *a, const void *b) { double x = fabs(*(const double*)a), y = fabs(*(const double*)b); return (x > y) - (x < y); }

extern void cblas_zdotu_sub(const CBLAS_INT N, const void *X, const CBLAS_INT incX, const void *Y, const CBLAS_INT incY, void *dotu);
extern void cblas_zdotu_sub_b(const CBLAS_INT N, const void *X, void *X_b, const CBLAS_INT incX, const void *Y, void *Y_b, const CBLAS_INT incY, void *dotu, void *dotu_b);

int main(void) {
    int i, j, idx, n_products;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofzx_(&diffblas_isize);
        set_isize1ofzy_(&diffblas_isize);
    }
    double h = 1.0e-7;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double vjp_fd, vjp_ad;

    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT incX = 1;
    CBLAS_INT incY = 1;

    double complex X[MAX_SIZE], X_b[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    double complex Y[MAX_SIZE], Y_b[MAX_SIZE], Y_orig[MAX_SIZE], Y_dir[MAX_SIZE];
    double complex dotu[1], dotu_b[1], dotu_orig[1], dotu_dir[1];
    double complex dotu_plus[1], dotu_minus[1], dotu_central_diff[1], dotu_b_orig[1];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { Y[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < 1; i++) { dotu[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); dotu_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }

    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));
    memcpy(Y_orig, Y, sizeof(Y[0])*(MAX_SIZE));
    memcpy(dotu_orig, dotu, sizeof(dotu[0])*(1));

    for (i = 0; i < 1; i++) { dotu_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); dotu_b_orig[i] = dotu_b[i]; }
    for (i = 0; i < MAX_SIZE; i++) X_b[i] = 0.0;
    for (i = 0; i < MAX_SIZE; i++) Y_b[i] = 0.0;

    cblas_zdotu_sub_b(n, X, X_b, incX, Y, Y_b, incY, dotu, dotu_b);

    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
    for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] + h * Y_dir[i];
    for (i = 0; i < 1; i++) dotu[i] = dotu_orig[i] + h * dotu_dir[i];
    cblas_zdotu_sub(n, X, incX, Y, incY, dotu);
    memcpy(dotu_plus, dotu, sizeof(dotu[0])*(1));

    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
    for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] - h * Y_dir[i];
    for (i = 0; i < 1; i++) dotu[i] = dotu_orig[i] - h * dotu_dir[i];
    cblas_zdotu_sub(n, X, incX, Y, incY, dotu);
    memcpy(dotu_minus, dotu, sizeof(dotu[0])*(1));

    vjp_fd = 0.0;
    {
        double temp_products[1];
        for (i = 0; i < 1; i++) temp_products[i] = creal(conj(dotu_b_orig[i]) * ((dotu_plus[i] - dotu_minus[i]) / (2.0*h)));
        qsort(temp_products, (size_t)1, sizeof(double), compare_abs_d);
        for (idx = 0; idx < 1; idx++) vjp_fd += temp_products[idx];
    }

    vjp_ad = 0.0;
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(X_dir[i]) * X_b[i]);
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(Y_dir[i]) * Y_b[i]);
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

