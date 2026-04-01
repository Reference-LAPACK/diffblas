/* Test program for cblas_dspr reverse mode (VJP verification, generic) */
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

static int compare_abs_d(const void *a, const void *b) { double x = fabs(*(const double*)a), y = fabs(*(const double*)b); return (x > y) - (x < y); }

extern void cblas_dspr(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_INT N, const double alpha, const double *X, const CBLAS_INT incX, double *Ap);
extern void cblas_dspr_b(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_INT N, const double alpha, double *alpha_b, const double *X, double *X_b, const CBLAS_INT incX, double *Ap, double *Ap_b);

int main(void) {
    int i, j, idx, n_products;
    double h = 1.0e-7;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT incX = 1;

    double alpha, alpha_b, alpha_orig, alpha_dir;
    double X[MAX_SIZE], X_b[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    double Ap[PACKED_SIZE], Ap_b[PACKED_SIZE], Ap_orig[PACKED_SIZE], Ap_dir[PACKED_SIZE];
    double Ap_plus[PACKED_SIZE], Ap_minus[PACKED_SIZE], Ap_central_diff[PACKED_SIZE], Ap_b_orig[PACKED_SIZE];

    srand(42);
    alpha = (rand()/(double)RAND_MAX)*2.0 - 1.0; alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < PACKED_SIZE; i++) { Ap[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; Ap_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }

    alpha_orig = alpha;
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));
    memcpy(Ap_orig, Ap, sizeof(Ap[0])*(PACKED_SIZE));

    for (i = 0; i < PACKED_SIZE; i++) { Ap_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; Ap_b_orig[i] = Ap_b[i]; }
    alpha_b = 0.0;
    for (i = 0; i < MAX_SIZE; i++) X_b[i] = 0.0;

    cblas_dspr_b(layout, uplo, n, alpha, &alpha_b, X, X_b, incX, Ap, Ap_b);

    alpha = alpha_orig + h * alpha_dir;
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
    for (i = 0; i < PACKED_SIZE; i++) Ap[i] = Ap_orig[i] + h * Ap_dir[i];
    cblas_dspr(layout, uplo, n, alpha, X, incX, Ap);
    memcpy(Ap_plus, Ap, sizeof(Ap[0])*(PACKED_SIZE));

    alpha = alpha_orig - h * alpha_dir;
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
    for (i = 0; i < PACKED_SIZE; i++) Ap[i] = Ap_orig[i] - h * Ap_dir[i];
    cblas_dspr(layout, uplo, n, alpha, X, incX, Ap);
    memcpy(Ap_minus, Ap, sizeof(Ap[0])*(PACKED_SIZE));

    vjp_fd = 0.0;
    {
        double temp_products[PACKED_SIZE];
        for (i = 0; i < PACKED_SIZE; i++) temp_products[i] = Ap_b_orig[i] * ((Ap_plus[i] - Ap_minus[i]) / (2.0*h));
        qsort(temp_products, (size_t)PACKED_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < PACKED_SIZE; idx++) vjp_fd += temp_products[idx];
    }

    vjp_ad = 0.0;
    vjp_ad += alpha_dir * alpha_b;
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = X_dir[i] * X_b[i];
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        double temp_products[PACKED_SIZE];
        for (i = 0; i < PACKED_SIZE; i++) temp_products[i] = Ap_dir[i] * Ap_b[i];
        qsort(temp_products, (size_t)PACKED_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < PACKED_SIZE; idx++) vjp_ad += temp_products[idx];
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

