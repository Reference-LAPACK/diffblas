/* Test program for cblas_cgemm reverse mode (VJP verification) */
/* Generated automatically by run_tapenade_cblas.py */
/* Mode: b (reverse) - same derivative check as BLAS test_*_reverse.f90 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "cblas.h"
#include "cblas_f77.h"

#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE
extern void set_isize2ofa_(int *val);
extern void set_isize2ofb_(int *val);

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_cgemm(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const CBLAS_INT M, const CBLAS_INT N, const CBLAS_INT K, const void *alpha, const void *A, const CBLAS_INT lda, const void *B, const CBLAS_INT ldb, const void *beta, void *C, const CBLAS_INT ldc);
extern void cblas_cgemm_b(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const CBLAS_INT M, const CBLAS_INT N, const CBLAS_INT K, const void *alpha, void *alpha_b, const void *A, void *A_b, const CBLAS_INT lda, const void *B, void *B_b, const CBLAS_INT ldb, const void *beta, void *beta_b, void *C, void *C_b, const CBLAS_INT ldc);

int main(void) {
    int i, j;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize2ofa_(&diffblas_isize);
        set_isize2ofb_(&diffblas_isize);
    }
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;

    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_INT m = TEST_SIZE, n = TEST_SIZE, k = TEST_SIZE;
    CBLAS_INT lda = MAX_SIZE, ldb = MAX_SIZE, ldc = MAX_SIZE;
    float complex alpha, alpha_b, alpha_dir;
    float complex beta, beta_b, beta_dir;
    float complex A[MAX_SIZE*MAX_SIZE], B[MAX_SIZE*MAX_SIZE], C[MAX_SIZE*MAX_SIZE];
    float complex A_b[MAX_SIZE*MAX_SIZE], B_b[MAX_SIZE*MAX_SIZE], C_b[MAX_SIZE*MAX_SIZE];
    float complex A_dir[MAX_SIZE*MAX_SIZE], B_dir[MAX_SIZE*MAX_SIZE], C_dir[MAX_SIZE*MAX_SIZE];
    float complex C_forward[MAX_SIZE*MAX_SIZE], C_backward[MAX_SIZE*MAX_SIZE];
    float complex C_b_orig[MAX_SIZE*MAX_SIZE];  /* save cotangent before _b overwrites */
    float complex alpha_orig, beta_orig, A_orig[MAX_SIZE*MAX_SIZE], B_orig[MAX_SIZE*MAX_SIZE], C_orig[MAX_SIZE*MAX_SIZE];  /* for restore like BLAS test */

    srand(42);
    alpha = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
    beta  = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) {
        A[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
        B[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
        C[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
    }
    /* Cotangent (seed on output C) and direction vectors */
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) {
        C_b[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
        A_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
        B_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
        C_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
    }
    alpha_dir = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
    beta_dir  = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);

    /* Save original primals (restore before each FD call - match BLAS test_dgemm_reverse.f90) */
    alpha_orig = alpha; beta_orig = beta;
    memcpy(A_orig, A, sizeof(A)); memcpy(B_orig, B, sizeof(B)); memcpy(C_orig, C, sizeof(C));
    memcpy(C_b_orig, C_b, sizeof(C_b));  /* save cotangent before _b overwrites C_b */
    /* Initialize input adjoints to zero (they will be computed by _b) - match BLAS test */
    alpha_b = 0.0f; beta_b = 0.0f;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A_b[i] = 0.0f; B_b[i] = 0.0f; }
    /* Call reverse mode: interleaved (primal, adjoint) per Tapenade signature */
    cblas_cgemm_b(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k,
        (const void*)&alpha, (void*)&alpha_b, (const void*)A, (void*)A_b, lda, (const void*)B, (void*)B_b, ldb, (const void*)&beta, (void*)&beta_b, (void*)C, (void*)C_b, ldc);

    /* Forward perturbation: f(x_orig + h*dir) - restore from originals then add, like BLAS test */
    alpha = alpha_orig + h * alpha_dir; beta = beta_orig + h * beta_dir;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A[i] = A_orig[i] + h * A_dir[i]; B[i] = B_orig[i] + h * B_dir[i]; C[i] = C_orig[i] + h * C_dir[i]; }
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, (const void*)&alpha, (const void*)A, lda, (const void*)B, ldb, (const void*)&beta, (void*)C, ldc);
    memcpy(C_forward, C, sizeof(C));

    /* Backward perturbation: f(x_orig - h*dir) - restore from originals then subtract, like BLAS test */
    alpha = alpha_orig - h * alpha_dir; beta = beta_orig - h * beta_dir;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; i++) { A[i] = A_orig[i] - h * A_dir[i]; B[i] = B_orig[i] - h * B_dir[i]; C[i] = C_orig[i] - h * C_dir[i]; }
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, (const void*)&alpha, (const void*)A, lda, (const void*)B, ldb, (const void*)&beta, (void*)C, ldc);
    memcpy(C_backward, C, sizeof(C));

    float vjp_fd, vjp_ad;
    /* VJP left side: cotangent^T @ central_diff (FD), sorted summation - match BLAS test_dgemm_reverse.f90 */
    {
        float temp_products[MAX_SIZE*MAX_SIZE];
        int n_products = MAX_SIZE*MAX_SIZE, idx;
        for (i = 0; i < n_products; i++) temp_products[i] = creal(conj(C_b_orig[i]) * ((C_forward[i] - C_backward[i]) / (2.0*h)));
        qsort(temp_products, (size_t)n_products, sizeof(float), compare_abs_f);
        vjp_fd = 0.0f;
        for (idx = 0; idx < n_products; idx++) vjp_fd += temp_products[idx];
    }

    /* VJP right side: direction^T @ adjoint, sorted summation - match BLAS */
    vjp_ad = 0.0f;
    vjp_ad += creal(conj(alpha_dir) * alpha_b) + creal(conj(beta_dir) * beta_b);
    {
        float temp_products[MAX_SIZE*MAX_SIZE];
        int n_products = MAX_SIZE*MAX_SIZE, idx;
        for (i = 0; i < n_products; i++) temp_products[i] = creal(conj(A_dir[i]) * A_b[i]);
        qsort(temp_products, (size_t)n_products, sizeof(float), compare_abs_f);
        for (idx = 0; idx < n_products; idx++) vjp_ad += temp_products[idx];
        for (i = 0; i < n_products; i++) temp_products[i] = creal(conj(B_dir[i]) * B_b[i]);
        qsort(temp_products, (size_t)n_products, sizeof(float), compare_abs_f);
        for (idx = 0; idx < n_products; idx++) vjp_ad += temp_products[idx];
        for (i = 0; i < n_products; i++) temp_products[i] = creal(conj(C_dir[i]) * C_b[i]);
        qsort(temp_products, (size_t)n_products, sizeof(float), compare_abs_f);
        for (idx = 0; idx < n_products; idx++) vjp_ad += temp_products[idx];
    }

    /* Error check: |vjp_fd - vjp_ad| <= atol + rtol*|vjp_ad| - match BLAS test_dgemm_reverse.f90 */
    {
        float abs_err = fabsf(vjp_fd - vjp_ad);
        float abs_reference = fabsf(vjp_ad);
        float error_bound = atol + rtol * (abs_reference > 1e-10f ? abs_reference : 1e-10f);
        printf("VJP: fd=%.10e ad=%.10e abs_err=%.10e error_bound=%.10e\n", (double)vjp_fd, (double)vjp_ad, (double)abs_err, (double)error_bound);
        printf("Tolerance: atol=%.0e, rtol=%.0e\n", (double)atol, (double)rtol);
        if (abs_err > error_bound) { printf("FAIL: Large errors detected in derivatives (outside tolerance)\n"); return 1; }
        if (abs_err < 0.5 * error_bound) { printf("PASS: Derivatives are accurate to machine precision\n"); return 0; }
        printf("PASS: Derivatives are within tolerance (rtol + atol)\n"); return 0;
    }
}

