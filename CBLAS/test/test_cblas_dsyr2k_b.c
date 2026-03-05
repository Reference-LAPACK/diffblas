/* Test program for cblas_dsyr2k reverse mode (VJP verification, generic) */
/* Generated automatically by run_tapenade_cblas.py - same methodology as BLAS test_*_reverse.f90 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cblas.h"
#include "cblas_f77.h"

#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE
extern void set_isize2ofa_(int *val);
extern void set_isize2ofb_(int *val);

static int compare_abs_d(const void *a, const void *b) { double x = fabs(*(const double*)a), y = fabs(*(const double*)b); return (x > y) - (x < y); }

extern void cblas_dsyr2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE Trans, const CBLAS_INT N, const CBLAS_INT K, const double alpha, const double *A, const CBLAS_INT lda, const double *B, const CBLAS_INT ldb, const double beta, double *C, const CBLAS_INT ldc);
extern void cblas_dsyr2k_b(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE Trans, const CBLAS_INT N, const CBLAS_INT K, const double alpha, double *alpha_b, const double *A, double *A_b, const CBLAS_INT lda, const double *B, double *B_b, const CBLAS_INT ldb, const double beta, double *beta_b, double *C, double *C_b, const CBLAS_INT ldc);

int main(void) {
    int i, j, idx, n_products;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize2ofa_(&diffblas_isize);
        set_isize2ofb_(&diffblas_isize);
    }
    double h = 1.0e-7;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT k = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT ldb = MAX_SIZE;
    CBLAS_INT ldc = MAX_SIZE;

    double alpha, alpha_b, alpha_orig, alpha_dir;
    double A[MAX_SIZE*MAX_SIZE], A_b[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    double B[MAX_SIZE*MAX_SIZE], B_b[MAX_SIZE*MAX_SIZE], B_orig[MAX_SIZE*MAX_SIZE], B_dir[MAX_SIZE*MAX_SIZE];
    double beta, beta_b, beta_orig, beta_dir;
    double C[MAX_SIZE*MAX_SIZE], C_b[MAX_SIZE*MAX_SIZE], C_orig[MAX_SIZE*MAX_SIZE], C_dir[MAX_SIZE*MAX_SIZE];
    double C_plus[MAX_SIZE*MAX_SIZE], C_minus[MAX_SIZE*MAX_SIZE], C_central_diff[MAX_SIZE*MAX_SIZE], C_b_orig[MAX_SIZE*MAX_SIZE];

    srand(42);
    alpha = (rand()/(double)RAND_MAX)*2.0 - 1.0; alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { B[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; B_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    beta = (rand()/(double)RAND_MAX)*2.0 - 1.0; beta_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { C[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; C_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }

    alpha_orig = alpha;
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
    memcpy(B_orig, B, sizeof(B[0])*(MAX_SIZE*MAX_SIZE));
    beta_orig = beta;
    memcpy(C_orig, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));

    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { C_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; C_b_orig[i] = C_b[i]; }
    alpha_b = 0.0;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A_b[i] = 0.0;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) B_b[i] = 0.0;
    beta_b = 0.0;

    cblas_dsyr2k_b(layout, uplo, trans, n, k, alpha, &alpha_b, A, A_b, lda, B, B_b, ldb, beta, &beta_b, C, C_b, ldc);

    alpha = alpha_orig + h * alpha_dir;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] + h * A_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) B[i] = B_orig[i] + h * B_dir[i];
    beta = beta_orig + h * beta_dir;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) C[i] = C_orig[i] + h * C_dir[i];
    cblas_dsyr2k(layout, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    memcpy(C_plus, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));

    alpha = alpha_orig - h * alpha_dir;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) B[i] = B_orig[i] - h * B_dir[i];
    beta = beta_orig - h * beta_dir;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) C[i] = C_orig[i] - h * C_dir[i];
    cblas_dsyr2k(layout, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    memcpy(C_minus, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));

    vjp_fd = 0.0;
    {
        double temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = C_b_orig[i] * ((C_plus[i] - C_minus[i]) / (2.0*h));
        qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_fd += temp_products[idx];
    }

    vjp_ad = 0.0;
    vjp_ad += alpha_dir * alpha_b;
    {
        double temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = A_dir[i] * A_b[i];
        qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        double temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = B_dir[i] * B_b[i];
        qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    vjp_ad += beta_dir * beta_b;
    {
        double temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = C_dir[i] * C_b[i];
        qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
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

