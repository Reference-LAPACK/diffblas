/* Test program for cblas_csyr2k reverse mode (VJP verification, generic) */
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

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_csyr2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE Trans, const CBLAS_INT N, const CBLAS_INT K, const void *alpha, const void *A, const CBLAS_INT lda, const void *B, const CBLAS_INT ldb, const void *beta, void *C, const CBLAS_INT ldc);
extern void cblas_csyr2k_b(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE Trans, const CBLAS_INT N, const CBLAS_INT K, const void *alpha, void *alpha_b, const void *A, void *A_b, const CBLAS_INT lda, const void *B, void *B_b, const CBLAS_INT ldb, const void *beta, void *beta_b, void *C, void *C_b, const CBLAS_INT ldc);

int main(void) {
    int i, j, idx, n_products;
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;
    float vjp_fd, vjp_ad;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    CBLAS_UPLO uplo = CblasUpper;
    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT k = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE;
    CBLAS_INT ldb = MAX_SIZE;
    CBLAS_INT ldc = MAX_SIZE;

    float complex alpha[MAX_SIZE], alpha_b[MAX_SIZE], alpha_orig[MAX_SIZE], alpha_dir[MAX_SIZE];
    float complex A[MAX_SIZE*MAX_SIZE], A_b[MAX_SIZE*MAX_SIZE], A_orig[MAX_SIZE*MAX_SIZE], A_dir[MAX_SIZE*MAX_SIZE];
    float complex B[MAX_SIZE*MAX_SIZE], B_b[MAX_SIZE*MAX_SIZE], B_orig[MAX_SIZE*MAX_SIZE], B_dir[MAX_SIZE*MAX_SIZE];
    float complex beta[MAX_SIZE], beta_b[MAX_SIZE], beta_orig[MAX_SIZE], beta_dir[MAX_SIZE];
    float complex C[MAX_SIZE*MAX_SIZE], C_b[MAX_SIZE*MAX_SIZE], C_orig[MAX_SIZE*MAX_SIZE], C_dir[MAX_SIZE*MAX_SIZE];
    float complex C_plus[MAX_SIZE*MAX_SIZE], C_minus[MAX_SIZE*MAX_SIZE], C_central_diff[MAX_SIZE*MAX_SIZE], C_b_orig[MAX_SIZE*MAX_SIZE];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) { alpha[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); alpha_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { B[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); B_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { beta[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); beta_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { C[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); C_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }

    memcpy(alpha_orig, alpha, sizeof(alpha[0])*(MAX_SIZE));
    memcpy(A_orig, A, sizeof(A[0])*(MAX_SIZE*MAX_SIZE));
    memcpy(B_orig, B, sizeof(B[0])*(MAX_SIZE*MAX_SIZE));
    memcpy(beta_orig, beta, sizeof(beta[0])*(MAX_SIZE));
    memcpy(C_orig, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));

    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { C_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); C_b_orig[i] = C_b[i]; }
    for (i = 0; i < MAX_SIZE; i++) alpha_b[i] = 0.0f;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A_b[i] = 0.0f;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) B_b[i] = 0.0f;
    for (i = 0; i < MAX_SIZE; i++) beta_b[i] = 0.0f;

    cblas_csyr2k_b(layout, uplo, trans, n, k, &alpha, alpha_b, A, A_b, lda, B, B_b, ldb, &beta, beta_b, C, C_b, ldc);

    for (i = 0; i < MAX_SIZE; i++) alpha[i] = alpha_orig[i] + h * alpha_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] + h * A_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) B[i] = B_orig[i] + h * B_dir[i];
    for (i = 0; i < MAX_SIZE; i++) beta[i] = beta_orig[i] + h * beta_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) C[i] = C_orig[i] + h * C_dir[i];
    cblas_csyr2k(layout, uplo, trans, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
    memcpy(C_plus, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));

    for (i = 0; i < MAX_SIZE; i++) alpha[i] = alpha_orig[i] - h * alpha_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) A[i] = A_orig[i] - h * A_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) B[i] = B_orig[i] - h * B_dir[i];
    for (i = 0; i < MAX_SIZE; i++) beta[i] = beta_orig[i] - h * beta_dir[i];
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) C[i] = C_orig[i] - h * C_dir[i];
    cblas_csyr2k(layout, uplo, trans, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
    memcpy(C_minus, C, sizeof(C[0])*(MAX_SIZE*MAX_SIZE));

    vjp_fd = 0.0f;
    {
        float temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = creal(conj(C_b_orig[i]) * ((C_plus[i] - C_minus[i]) / (2.0*h)));
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
        float temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = creal(conj(A_dir[i]) * A_b[i]);
        qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        float temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = creal(conj(B_dir[i]) * B_b[i]);
        qsort(temp_products, (size_t)MAX_SIZE*MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE*MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        float temp_products[MAX_SIZE];
        for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(beta_dir[i]) * beta_b[i]);
        qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
        for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
    }
    {
        float temp_products[MAX_SIZE*MAX_SIZE];
        for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) temp_products[i] = creal(conj(C_dir[i]) * C_b[i]);
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

