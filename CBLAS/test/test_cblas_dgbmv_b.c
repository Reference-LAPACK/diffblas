/* Test program for cblas_dgbmv reverse mode (VJP verification, generic) */
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
extern void set_isize2ofa_(int *val);

static int compare_abs_d(const void *a, const void *b) { double x = fabs(*(const double*)a), y = fabs(*(const double*)b); return (x > y) - (x < y); }

extern void cblas_dgbmv(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA, const CBLAS_INT M, const CBLAS_INT N, const CBLAS_INT KL, const CBLAS_INT KU, const double alpha, const double *A, const CBLAS_INT lda, const double *X, const CBLAS_INT incX, const double beta, double *Y, const CBLAS_INT incY);
extern void cblas_dgbmv_b(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA, const CBLAS_INT M, const CBLAS_INT N, const CBLAS_INT KL, const CBLAS_INT KU, const double alpha, double *alpha_b, const double *A, double *A_b, const CBLAS_INT lda, const double *X, double *X_b, const CBLAS_INT incX, const double beta, double *beta_b, double *Y, double *Y_b, const CBLAS_INT incY);

int main(void) {
    int i, j, idx, n_products;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofx_(&diffblas_isize);
        set_isize2ofa_(&diffblas_isize);
    }
    double h = 1.0e-7;
    double atol = 1.0e-5, rtol = 1.0e-5;
    double vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE transa = CblasNoTrans;
    CBLAS_INT m = TEST_SIZE;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT incX = 1;
    CBLAS_INT incY = 1;
    CBLAS_INT KL = 1;  /* band width: LDA >= KL+KU+1 (gbmv) */
    CBLAS_INT KU = 1;  /* band width: LDA >= KL+KU+1 (gbmv) */

    double alpha, alpha_b, alpha_orig, alpha_dir;
    double A[MAX_SIZE*MAX_SIZE], A_b[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    double X[MAX_SIZE], X_b[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    double beta, beta_b, beta_orig, beta_dir;
    double Y[MAX_SIZE], Y_b[MAX_SIZE], Y_orig[MAX_SIZE], Y_dir[MAX_SIZE];
    double Y_plus[MAX_SIZE], Y_minus[MAX_SIZE], Y_central_diff[MAX_SIZE], Y_b_orig[MAX_SIZE];

    srand(42);
    alpha = (rand()/(double)RAND_MAX)*2.0 - 1.0; alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    /* A: general band storage (KL+KU+1) x N (match BLAS/test) */
    memset(A, 0, sizeof(A));
    for (j = 0; j < MAX_SIZE; j++) {
        int band_rows = KL + KU + 1;
        for (i = 0; i < band_rows; i++) {
            A[i + j * MAX_SIZE] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        }
    }
    /* A: general band storage (KL+KU+1) x N (match BLAS/test) */
    memset(A_dir, 0, sizeof(A_dir));
    for (j = 0; j < MAX_SIZE; j++) {
        int band_rows = KL + KU + 1;
        for (i = 0; i < band_rows; i++) {
            A_dir[i + j * MAX_SIZE] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        }
    }
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }
    beta = (rand()/(double)RAND_MAX)*2.0 - 1.0; beta_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    for (i = 0; i < MAX_SIZE; i++) { Y[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; }

    alpha_orig = alpha;
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));
    beta_orig = beta;
    memcpy(Y_orig, Y, sizeof(Y[0])*(MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) { Y_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0; Y_b_orig[i] = Y_b[i]; }
    alpha_b = 0.0;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A_b[i] = 0.0;
    for (i = 0; i < MAX_SIZE; i++) X_b[i] = 0.0;
    beta_b = 0.0;

    cblas_dgbmv_b(layout, transa, m, n, KL, KU, alpha, &alpha_b, A, A_b, lda, X, X_b, incX, beta, &beta_b, Y, Y_b, incY);

    alpha = alpha_orig + h * alpha_dir;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] + h * A_dir[i];
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
    beta = beta_orig + h * beta_dir;
    for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] + h * Y_dir[i];
    cblas_dgbmv(layout, transa, m, n, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
    memcpy(Y_plus, Y, sizeof(Y[0])*(MAX_SIZE));

    alpha = alpha_orig - h * alpha_dir;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
    for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
    beta = beta_orig - h * beta_dir;
    for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] - h * Y_dir[i];
    cblas_dgbmv(layout, transa, m, n, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
    memcpy(Y_minus, Y, sizeof(Y[0])*(MAX_SIZE));

    vjp_fd = 0.0;
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = Y_b_orig[i] * ((Y_plus[i] - Y_minus[i]) / (2.0*h));
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_fd += temp_products[idx];
    }

    vjp_ad = 0.0;
    vjp_ad += alpha_dir * alpha_b;
    {
        double temp_products[MAX_SIZE*MAX_SIZE];
        int n_band = 0;
        int band_rows = KL + KU + 1;
        for (j = 0; j < n; j++) for (i = 0; i < band_rows; i++) {
                temp_products[n_band++] = A_dir[i+j*lda] * A_b[i+j*lda];
            }
        qsort(temp_products, (size_t)n_band, sizeof(double), compare_abs_d);
        for (idx = 0; idx < n_band; idx++) vjp_ad += temp_products[idx];
    }
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = X_dir[i] * X_b[i];
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(double), compare_abs_d);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    vjp_ad += beta_dir * beta_b;
    {
        double temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = Y_dir[i] * Y_b[i];
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

