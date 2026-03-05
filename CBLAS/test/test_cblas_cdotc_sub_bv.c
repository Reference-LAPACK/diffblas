/* Test program for cblas_cdotc_sub vector reverse mode (VJP verification, generic, loop over directions) */
/* Generated automatically by run_tapenade_cblas.py */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "cblas.h"
#include "cblas_f77.h"
#include "cblas_bv.h"

#ifndef NBDirsMax
#define NBDirsMax 4
#endif
#define TEST_SIZE 4
#define MAX_SIZE TEST_SIZE
extern void set_isize1ofcx_(int *val);
extern void set_isize1ofcy_(int *val);

static int compare_abs_f(const void *a, const void *b) { float x = fabsf(*(const float*)a), y = fabsf(*(const float*)b); return (x > y) - (x < y); }

extern void cblas_cdotc_sub(const CBLAS_INT N, const void *X, const CBLAS_INT incX, const void *Y, const CBLAS_INT incY, void *dotc);
/* cblas_*_bv from cblas_bv.h */

int main(void) {
    int i, j, idx, idir, nbdirs = NBDirsMax, n_products;
    {
        int diffblas_isize = MAX_SIZE;
        set_isize1ofcx_(&diffblas_isize);
        set_isize1ofcy_(&diffblas_isize);
    }
    int has_large_errors = 0;
    float h = 1.0e-3f;
    float atol = 1.0e-2f, rtol = 1.0e-2f;
    float max_error = 0.0f;
    float vjp_fd, vjp_ad;

    CBLAS_INT n = TEST_SIZE;
    CBLAS_INT incX = 1;
    CBLAS_INT incY = 1;

    float complex X[MAX_SIZE], X_orig[MAX_SIZE], X_dir[MAX_SIZE];
    float complex X_b[MAX_SIZE][NBDirsMax], X_b_orig[MAX_SIZE][NBDirsMax];
    float complex Y[MAX_SIZE], Y_orig[MAX_SIZE], Y_dir[MAX_SIZE];
    float complex Y_b[MAX_SIZE][NBDirsMax], Y_b_orig[MAX_SIZE][NBDirsMax];
    float complex dotc[1], dotc_orig[1], dotc_dir[1];
    float complex dotc_b[1][NBDirsMax], dotc_b_orig[1][NBDirsMax];
    float complex dotc_plus[1], dotc_minus[1];

    srand(42);
    for (i = 0; i < MAX_SIZE; i++) { X[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < MAX_SIZE; i++) { Y[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }
    for (i = 0; i < 1; i++) { dotc[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); dotc_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); }

    memcpy(X_orig, X, sizeof(X[0])*(MAX_SIZE));
    memcpy(Y_orig, Y, sizeof(Y[0])*(MAX_SIZE));
    memcpy(dotc_orig, dotc, sizeof(dotc[0])*(1));

    for (i = 0; i < 1; i++) for (j = 0; j < NBDirsMax; j++) { dotc_b[i][j] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0); dotc_b_orig[i][j] = dotc_b[i][j]; }
    for (i = 0; i < MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) X_b[i][j] = 0.0f;
    for (i = 0; i < MAX_SIZE; i++) for (j = 0; j < NBDirsMax; j++) Y_b[i][j] = 0.0f;

    cblas_cdotc_sub_bv(n, X, (void*)X_b, incX, Y, (void*)Y_b, incY, dotc, (void*)dotc_b, nbdirs);

    for (idir = 0; idir < nbdirs; idir++) {
        /* Restore primals for this direction */
        memcpy(X, X_orig, sizeof(X[0])*(MAX_SIZE));
        memcpy(Y, Y_orig, sizeof(Y[0])*(MAX_SIZE));
        memcpy(dotc, dotc_orig, sizeof(dotc[0])*(1));
        /* Random direction for this idir */
        for (i = 0; i < MAX_SIZE; i++) X_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
        for (i = 0; i < MAX_SIZE; i++) Y_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
        for (i = 0; i < 1; i++) dotc_dir[i] = (rand()/(double)RAND_MAX)*2.0 - 1.0 + I*((rand()/(double)RAND_MAX)*2.0 - 1.0);
        /* Forward */
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] + h * X_dir[i];
        for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] + h * Y_dir[i];
        for (i = 0; i < 1; i++) dotc[i] = dotc_orig[i] + h * dotc_dir[i];
        cblas_cdotc_sub(n, X, incX, Y, incY, dotc);
        memcpy(dotc_plus, dotc, sizeof(dotc[0])*(1));
        /* Backward */
        for (i = 0; i < MAX_SIZE; i++) X[i] = X_orig[i] - h * X_dir[i];
        for (i = 0; i < MAX_SIZE; i++) Y[i] = Y_orig[i] - h * Y_dir[i];
        for (i = 0; i < 1; i++) dotc[i] = dotc_orig[i] - h * dotc_dir[i];
        cblas_cdotc_sub(n, X, incX, Y, incY, dotc);
        memcpy(dotc_minus, dotc, sizeof(dotc[0])*(1));

        vjp_fd = 0.0f;
        {
            float temp_products[1];
            for (i = 0; i < 1; i++) temp_products[i] = creal(conj(dotc_b_orig[i][idir]) * ((dotc_plus[i] - dotc_minus[i]) / (2.0*h)));
            qsort(temp_products, (size_t)1, sizeof(float), compare_abs_f);
            for (idx = 0; idx < 1; idx++) vjp_fd += temp_products[idx];
        }
        vjp_ad = 0.0f;
        {
            float temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(X_dir[i]) * X_b[i][idir]);
            qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }
        {
            float temp_products[MAX_SIZE];
            for (i = 0; i < MAX_SIZE; i++) temp_products[i] = creal(conj(Y_dir[i]) * Y_b[i][idir]);
            qsort(temp_products, (size_t)MAX_SIZE, sizeof(float), compare_abs_f);
            for (idx = 0; idx < MAX_SIZE; idx++) vjp_ad += temp_products[idx];
        }

        {
            float abs_err = fabsf(vjp_fd - vjp_ad);
            float abs_reference = fabsf(vjp_ad);
            float error_bound = atol + rtol * (abs_reference > 1e-10f ? abs_reference : 1e-10f);
            if (abs_err > error_bound) has_large_errors = 1;
            { float r = abs_err / error_bound; if (r > max_error) max_error = r; }
        }
    }

    printf("Maximum error ratio: %.6e\n", (double)max_error);
    if (has_large_errors) { printf("FAIL: Large errors in derivatives\n"); return 1; }
    if (max_error < 0.5f) { printf("PASS: Derivatives accurate to machine precision\n"); return 0; }
    if (max_error < 1.0f) { printf("PASS: Derivatives reasonably accurate\n"); return 0; }
    printf("WARNING: Derivatives may have significant errors\n"); return 0;
}

