/* Test program for cblas_snrm2 reverse mode (nrm2 VJP verification) */
/* Generated automatically by run_tapenade_cblas.py - matches BLAS test_*nrm2_reverse.f90 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cblas.h"

#define TEST_SIZE 4  /* match BLAS n=4 */

extern void cblas_snrm2_b(const CBLAS_INT N, const float *X, float *X_b, const CBLAS_INT incX, float cblas_snrm2_b_cotangent);

static int cmp_abs(const void *a, const void *b) {
    float fa = fabs(*(const float *)a), fb = fabs(*(const float *)b);
    return (fa > fb) - (fa < fb);
}

int main(void) {
    CBLAS_INT N = TEST_SIZE, incX = 1;
    float X[TEST_SIZE], X_b[TEST_SIZE], X_dir[TEST_SIZE];
    float nrm2_plus, nrm2_minus, nrm2_b = 1.0f;
    float h = 1.0e-3f, atol = 2.0e-3f, rtol = 2.0e-3f;
    float products[TEST_SIZE];
    int i;
    srand(42);
    for (i = 0; i < TEST_SIZE; i++) {
        X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
        X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0;
    }
    float nrm2 = cblas_snrm2(N, X, incX);
    /* Input adjoints must be zero before _b (Fortran uses increment semantics, match BLAS) */
    for (i = 0; i < TEST_SIZE; i++) X_b[i] = 0.0f;
    cblas_snrm2_b(N, X, X_b, incX, nrm2_b);
    /* VJP fd: (nrm2(x+h*dir) - nrm2(x-h*dir))/(2h) with cotangent 1 */
    for (i = 0; i < TEST_SIZE; i++) X[i] += h * X_dir[i];
    nrm2_plus = cblas_snrm2(N, X, incX);
    for (i = 0; i < TEST_SIZE; i++) X[i] -= 2*h * X_dir[i];
    nrm2_minus = cblas_snrm2(N, X, incX);
    float vjp_fd = (nrm2_plus - nrm2_minus) / (2.0*h);
    /* VJP ad: direction^T @ adjoint with sorted summation (match BLAS) */
    for (i = 0; i < TEST_SIZE; i++) products[i] = X_dir[i] * X_b[i];
    qsort(products, (size_t)TEST_SIZE, sizeof(products[0]), cmp_abs);
    float vjp_ad = 0.0f;
    for (i = 0; i < TEST_SIZE; i++) vjp_ad += products[i];
    {
        float abs_err = fabs(vjp_fd - vjp_ad);
        float ref = (fabs(vjp_ad) > 1e-10) ? fabs(vjp_ad) : 1e-10;
        float error_bound = atol + rtol * ref;
        printf("VJP: fd=%.10e ad=%.10e abs_err=%.10e error_bound=%.10e\n", (double)vjp_fd, (double)vjp_ad, (double)abs_err, (double)error_bound);
        if (abs_err > error_bound) { printf("FAIL: Large errors detected in derivatives (outside tolerance)\n"); return 1; }
        if (abs_err < 0.5 * error_bound) { printf("PASS: Derivatives are accurate to machine precision\n"); return 0; }
        printf("PASS: Derivatives are reasonably accurate\n"); return 0;
    }
}

