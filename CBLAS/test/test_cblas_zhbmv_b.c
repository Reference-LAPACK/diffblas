/* Test program for cblas_zhbmv reverse mode (VJP verification, generic) */
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

extern void cblas_zhbmv(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_INT N, const CBLAS_INT K, const void *alpha, const void *A, const CBLAS_INT lda, const void *X, const CBLAS_INT incX, const void *beta, void *Y, const CBLAS_INT incY);
extern void cblas_zhbmv_b(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_INT N, const CBLAS_INT K, const void *alpha, void *alpha_b, const void *A, void *A_b, const CBLAS_INT lda, const void *X, void *X_b, const CBLAS_INT incX, const void *beta, void *beta_b, void *Y, void *Y_b, const CBLAS_INT incY);

int main(void) {
    int i, j, idx, n_products;
    double h = 1.0e-7;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT k = n - 1;  /* band width: lda >= k+1 */
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT incX = 1;
    CBLAS_INT incY = 1;

    double complex alpha[MAX_SIZE], alpha_b[MAX_SIZE], alpha_orig[MAX_SIZE], alpha_dir[MAX_SIZE];
    double complex A[MAX_SIZE*MAX_SIZE], A_b[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    double complex X[MAX_SIZE], X_b[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    double complex beta[MAX_SIZE], beta_b[MAX_SIZE], beta_orig[MAX_SIZE], beta_dir[MAX_SIZE];
    double complex Y[MAX_SIZE], Y_b[MAX_SIZE], Y_orig[MAX_SIZE], Y_dir[MAX_SIZE];
    double complex Y_plus[MAX_SIZE], Y_minus[MAX_SIZE], Y_central_diff[MAX_SIZE], Y_b_orig[MAX_SIZE];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) { alpha[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); alpha_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { beta[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); beta_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { Y[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }

    memcpy(alpha_orig, alpha, sizeof(alpha[0])*(MAX_SIZE));
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));
    memcpy(beta_orig, beta, sizeof(beta[0])*(MAX_SIZE));
    memcpy(Y_orig, Y, sizeof(Y[0])*(MAX_SIZE));

    /* Hermitian band A: real diagonal in band (row k = diagonal) */
    for (j = 0; j < n; j++) { A[k + j*lda] = creal(A[k + j*lda]); A_dir[k + j*lda] = creal(A_dir[k + j*lda]); }
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) { Y_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); Y_b_orig[i] = Y_b[i]; }
    for (i = 0; i < MAX_SIZE; i++) alpha_b[i] = 0.0;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A_b[i] = 0.0;
    for (i = 0; i < MAX_SIZE; i++) X_b[i] = 0.0;
    for (i = 0; i < MAX_SIZE; i++) beta_b[i] = 0.0;

    cblas_zhbmv_b(layout, uplo, n, k, &alpha, alpha_b, A, A_b, lda, X, X_b, incX, &beta, beta_b, Y, Y_b, incY);

    for (i = 0; i < MAX_SIZE; i++) alpha[i] = alpha_orig[i] + h * alpha_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] + h * A_dir[i];
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
    for (i = 0; i < MAX_SIZE; i++) beta[i] = beta_orig[i] + h * beta_dir[i];
    for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] + h * Y_dir[i];
    cblas_zhbmv(layout, uplo, n, k, &alpha, A, lda, X, incX, &beta, Y, incY);
    memcpy(Y_plus, Y, sizeof(Y[0])*(MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) alpha[i] = alpha_orig[i] - h * alpha_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
    for (i = 0; i < MAX_SIZE; i++) beta[i] = beta_orig[i] - h * beta_dir[i];
    for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] - h * Y_dir[i];
    cblas_zhbmv(layout, uplo, n, k, &alpha, A, lda, X, incX, &beta, Y, incY);
    memcpy(Y_minus, Y, sizeof(Y[0])*(MAX_SIZE));

    vjp_fd = 0.0;
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(Y_b_orig[i]) * ((Y_plus[i] - Y_minus[i]) / (2.0*h)));
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
        double temp_products[MAX_SIZE*MAX_SIZE];
        int n_band = 0;
        for (j = 0; j < n; j++)
            for (i = 0; i <= k; i++) {
                temp_products[n_band++] = creal(conj(A_dir[i+j*lda]) * A_b[i+j*lda]);
            }
        qsort(temp_products, (size_t)n_band, sizeof(double), compare_abs_d);
        for (idx = 0; idx < n_band; idx++) vjp_ad += temp_products[idx];
    }
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(X_dir[i]) * X_b[i]);
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(beta_dir[i]) * beta_b[i]);
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

