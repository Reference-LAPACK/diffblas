# 0 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/src/cblas_dger.c"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 0 "<command-line>" 2
# 1 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/src/cblas_dger.c"
# 10 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/src/cblas_dger.c"
# 1 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h" 1


# 1 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stddef.h" 1 3 4
# 145 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stddef.h" 3 4

# 145 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stddef.h" 3 4
typedef long int ptrdiff_t;
# 214 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stddef.h" 3 4
typedef long unsigned int size_t;
# 329 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stddef.h" 3 4
typedef int wchar_t;
# 425 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stddef.h" 3 4
typedef struct {
  long long __max_align_ll __attribute__((__aligned__(__alignof__(long long))));
  long double __max_align_ld __attribute__((__aligned__(__alignof__(long double))));
# 436 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stddef.h" 3 4
} max_align_t;
# 4 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h" 2
# 1 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stdint.h" 1 3 4
# 9 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stdint.h" 3 4
# 1 "/usr/include/stdint.h" 1 3 4
# 26 "/usr/include/stdint.h" 3 4
# 1 "/usr/include/bits/libc-header-start.h" 1 3 4
# 33 "/usr/include/bits/libc-header-start.h" 3 4
# 1 "/usr/include/features.h" 1 3 4
# 438 "/usr/include/features.h" 3 4
# 1 "/usr/include/sys/cdefs.h" 1 3 4
# 501 "/usr/include/sys/cdefs.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 502 "/usr/include/sys/cdefs.h" 2 3 4
# 1 "/usr/include/bits/long-double.h" 1 3 4
# 503 "/usr/include/sys/cdefs.h" 2 3 4
# 439 "/usr/include/features.h" 2 3 4
# 462 "/usr/include/features.h" 3 4
# 1 "/usr/include/gnu/stubs.h" 1 3 4
# 10 "/usr/include/gnu/stubs.h" 3 4
# 1 "/usr/include/gnu/stubs-64.h" 1 3 4
# 11 "/usr/include/gnu/stubs.h" 2 3 4
# 463 "/usr/include/features.h" 2 3 4
# 34 "/usr/include/bits/libc-header-start.h" 2 3 4
# 27 "/usr/include/stdint.h" 2 3 4
# 1 "/usr/include/bits/types.h" 1 3 4
# 27 "/usr/include/bits/types.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 28 "/usr/include/bits/types.h" 2 3 4


typedef unsigned char __u_char;
typedef unsigned short int __u_short;
typedef unsigned int __u_int;
typedef unsigned long int __u_long;


typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef signed short int __int16_t;
typedef unsigned short int __uint16_t;
typedef signed int __int32_t;
typedef unsigned int __uint32_t;

typedef signed long int __int64_t;
typedef unsigned long int __uint64_t;






typedef __int8_t __int_least8_t;
typedef __uint8_t __uint_least8_t;
typedef __int16_t __int_least16_t;
typedef __uint16_t __uint_least16_t;
typedef __int32_t __int_least32_t;
typedef __uint32_t __uint_least32_t;
typedef __int64_t __int_least64_t;
typedef __uint64_t __uint_least64_t;



typedef long int __quad_t;
typedef unsigned long int __u_quad_t;







typedef long int __intmax_t;
typedef unsigned long int __uintmax_t;
# 140 "/usr/include/bits/types.h" 3 4
# 1 "/usr/include/bits/typesizes.h" 1 3 4
# 141 "/usr/include/bits/types.h" 2 3 4


typedef unsigned long int __dev_t;
typedef unsigned int __uid_t;
typedef unsigned int __gid_t;
typedef unsigned long int __ino_t;
typedef unsigned long int __ino64_t;
typedef unsigned int __mode_t;
typedef unsigned long int __nlink_t;
typedef long int __off_t;
typedef long int __off64_t;
typedef int __pid_t;
typedef struct { int __val[2]; } __fsid_t;
typedef long int __clock_t;
typedef unsigned long int __rlim_t;
typedef unsigned long int __rlim64_t;
typedef unsigned int __id_t;
typedef long int __time_t;
typedef unsigned int __useconds_t;
typedef long int __suseconds_t;

typedef int __daddr_t;
typedef int __key_t;


typedef int __clockid_t;


typedef void * __timer_t;


typedef long int __blksize_t;




typedef long int __blkcnt_t;
typedef long int __blkcnt64_t;


typedef unsigned long int __fsblkcnt_t;
typedef unsigned long int __fsblkcnt64_t;


typedef unsigned long int __fsfilcnt_t;
typedef unsigned long int __fsfilcnt64_t;


typedef long int __fsword_t;

typedef long int __ssize_t;


typedef long int __syscall_slong_t;

typedef unsigned long int __syscall_ulong_t;



typedef __off64_t __loff_t;
typedef char *__caddr_t;


typedef long int __intptr_t;


typedef unsigned int __socklen_t;




typedef int __sig_atomic_t;
# 28 "/usr/include/stdint.h" 2 3 4
# 1 "/usr/include/bits/wchar.h" 1 3 4
# 29 "/usr/include/stdint.h" 2 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 30 "/usr/include/stdint.h" 2 3 4




# 1 "/usr/include/bits/stdint-intn.h" 1 3 4
# 24 "/usr/include/bits/stdint-intn.h" 3 4
typedef __int8_t int8_t;
typedef __int16_t int16_t;
typedef __int32_t int32_t;
typedef __int64_t int64_t;
# 35 "/usr/include/stdint.h" 2 3 4


# 1 "/usr/include/bits/stdint-uintn.h" 1 3 4
# 24 "/usr/include/bits/stdint-uintn.h" 3 4
typedef __uint8_t uint8_t;
typedef __uint16_t uint16_t;
typedef __uint32_t uint32_t;
typedef __uint64_t uint64_t;
# 38 "/usr/include/stdint.h" 2 3 4





typedef __int_least8_t int_least8_t;
typedef __int_least16_t int_least16_t;
typedef __int_least32_t int_least32_t;
typedef __int_least64_t int_least64_t;


typedef __uint_least8_t uint_least8_t;
typedef __uint_least16_t uint_least16_t;
typedef __uint_least32_t uint_least32_t;
typedef __uint_least64_t uint_least64_t;





typedef signed char int_fast8_t;

typedef long int int_fast16_t;
typedef long int int_fast32_t;
typedef long int int_fast64_t;
# 71 "/usr/include/stdint.h" 3 4
typedef unsigned char uint_fast8_t;

typedef unsigned long int uint_fast16_t;
typedef unsigned long int uint_fast32_t;
typedef unsigned long int uint_fast64_t;
# 87 "/usr/include/stdint.h" 3 4
typedef long int intptr_t;


typedef unsigned long int uintptr_t;
# 101 "/usr/include/stdint.h" 3 4
typedef __intmax_t intmax_t;
typedef __uintmax_t uintmax_t;
# 10 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stdint.h" 2 3 4
# 5 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h" 2
# 1 "/usr/include/inttypes.h" 1 3 4
# 34 "/usr/include/inttypes.h" 3 4
typedef int __gwchar_t;
# 266 "/usr/include/inttypes.h" 3 4





typedef struct
  {
    long int quot;
    long int rem;
  } imaxdiv_t;
# 290 "/usr/include/inttypes.h" 3 4
extern intmax_t imaxabs (intmax_t __n) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern imaxdiv_t imaxdiv (intmax_t __numer, intmax_t __denom)
      __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern intmax_t strtoimax (const char *__restrict __nptr,
      char **__restrict __endptr, int __base) __attribute__ ((__nothrow__ , __leaf__));


extern uintmax_t strtoumax (const char *__restrict __nptr,
       char ** __restrict __endptr, int __base) __attribute__ ((__nothrow__ , __leaf__));


extern intmax_t wcstoimax (const __gwchar_t *__restrict __nptr,
      __gwchar_t **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__));


extern uintmax_t wcstoumax (const __gwchar_t *__restrict __nptr,
       __gwchar_t ** __restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__));
# 432 "/usr/include/inttypes.h" 3 4

# 6 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h" 2
# 39 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h"

# 39 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h"
typedef enum CBLAS_LAYOUT {CblasRowMajor=101, CblasColMajor=102} CBLAS_LAYOUT;
typedef enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113} CBLAS_TRANSPOSE;
typedef enum CBLAS_UPLO {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
typedef enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132} CBLAS_DIAG;
typedef enum CBLAS_SIDE {CblasLeft=141, CblasRight=142} CBLAS_SIDE;



# 1 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas_mangling.h" 1
# 48 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h" 2
# 68 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h"
double cblas_dcabs1(const void *z);
float cblas_scabs1(const void *c);

float cblas_sdsdot(const int32_t N, const float alpha, const float *X,
                    const int32_t incX, const float *Y, const int32_t incY);
double cblas_dsdot(const int32_t N, const float *X, const int32_t incX, const float *Y,
                   const int32_t incY);
float cblas_sdot(const int32_t N, const float *X, const int32_t incX,
                  const float *Y, const int32_t incY);
double cblas_ddot(const int32_t N, const double *X, const int32_t incX,
                  const double *Y, const int32_t incY);




void cblas_cdotu_sub(const int32_t N, const void *X, const int32_t incX,
                       const void *Y, const int32_t incY, void *dotu);
void cblas_cdotc_sub(const int32_t N, const void *X, const int32_t incX,
                       const void *Y, const int32_t incY, void *dotc);

void cblas_zdotu_sub(const int32_t N, const void *X, const int32_t incX,
                       const void *Y, const int32_t incY, void *dotu);
void cblas_zdotc_sub(const int32_t N, const void *X, const int32_t incX,
                       const void *Y, const int32_t incY, void *dotc);





float cblas_snrm2(const int32_t N, const float *X, const int32_t incX);
float cblas_sasum(const int32_t N, const float *X, const int32_t incX);

double cblas_dnrm2(const int32_t N, const double *X, const int32_t incX);
double cblas_dasum(const int32_t N, const double *X, const int32_t incX);

float cblas_scnrm2(const int32_t N, const void *X, const int32_t incX);
float cblas_scasum(const int32_t N, const void *X, const int32_t incX);

double cblas_dznrm2(const int32_t N, const void *X, const int32_t incX);
double cblas_dzasum(const int32_t N, const void *X, const int32_t incX);





size_t cblas_isamax(const int32_t N, const float *X, const int32_t incX);
size_t cblas_idamax(const int32_t N, const double *X, const int32_t incX);
size_t cblas_icamax(const int32_t N, const void *X, const int32_t incX);
size_t cblas_izamax(const int32_t N, const void *X, const int32_t incX);
# 127 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h"
void cblas_sswap(const int32_t N, float *X, const int32_t incX,
                 float *Y, const int32_t incY);
void cblas_scopy(const int32_t N, const float *X, const int32_t incX,
                 float *Y, const int32_t incY);
void cblas_saxpy(const int32_t N, const float alpha, const float *X,
                 const int32_t incX, float *Y, const int32_t incY);

void cblas_dswap(const int32_t N, double *X, const int32_t incX,
                 double *Y, const int32_t incY);
void cblas_dcopy(const int32_t N, const double *X, const int32_t incX,
                 double *Y, const int32_t incY);
void cblas_daxpy(const int32_t N, const double alpha, const double *X,
                 const int32_t incX, double *Y, const int32_t incY);

void cblas_cswap(const int32_t N, void *X, const int32_t incX,
                 void *Y, const int32_t incY);
void cblas_ccopy(const int32_t N, const void *X, const int32_t incX,
                 void *Y, const int32_t incY);
void cblas_caxpy(const int32_t N, const void *alpha, const void *X,
                 const int32_t incX, void *Y, const int32_t incY);

void cblas_zswap(const int32_t N, void *X, const int32_t incX,
                 void *Y, const int32_t incY);
void cblas_zcopy(const int32_t N, const void *X, const int32_t incX,
                 void *Y, const int32_t incY);
void cblas_zaxpy(const int32_t N, const void *alpha, const void *X,
                 const int32_t incX, void *Y, const int32_t incY);





void cblas_srotmg(float *d1, float *d2, float *b1, const float b2, float *P);
void cblas_srotm(const int32_t N, float *X, const int32_t incX,
                 float *Y, const int32_t incY, const float *P);
void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double *P);
void cblas_drotm(const int32_t N, double *X, const int32_t incX,
                 double *Y, const int32_t incY, const double *P);






void cblas_sscal(const int32_t N, const float alpha, float *X, const int32_t incX);
void cblas_dscal(const int32_t N, const double alpha, double *X, const int32_t incX);
void cblas_cscal(const int32_t N, const void *alpha, void *X, const int32_t incX);
void cblas_zscal(const int32_t N, const void *alpha, void *X, const int32_t incX);
void cblas_csscal(const int32_t N, const float alpha, void *X, const int32_t incX);
void cblas_zdscal(const int32_t N, const double alpha, void *X, const int32_t incX);

void cblas_srotg(float *a, float *b, float *c, float *s);
void cblas_drotg(double *a, double *b, double *c, double *s);
void cblas_crotg(void *a, void *b, float *c, void *s);
void cblas_zrotg(void *a, void *b, double *c, void *s);

void cblas_srot(const int32_t N, float *X, const int32_t incX,
                float *Y, const int32_t incY, const float c, const float s);
void cblas_drot(const int32_t N, double *X, const int32_t incX,
                double *Y, const int32_t incY, const double c, const double s);
void cblas_csrot(const int32_t N, void *X, const int32_t incX,
                 void *Y, const int32_t incY, const float c, const float s);
void cblas_zdrot(const int32_t N, void *X, const int32_t incX,
                 void *Y, const int32_t incY, const double c, const double s);
# 201 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h"
void cblas_sgemv(const CBLAS_LAYOUT layout,
                 const CBLAS_TRANSPOSE TransA, const int32_t M, const int32_t N,
                 const float alpha, const float *A, const int32_t lda,
                 const float *X, const int32_t incX, const float beta,
                 float *Y, const int32_t incY);
void cblas_sgbmv(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int32_t M, const int32_t N,
                 const int32_t KL, const int32_t KU, const float alpha,
                 const float *A, const int32_t lda, const float *X,
                 const int32_t incX, const float beta, float *Y, const int32_t incY);
void cblas_strmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const float *A, const int32_t lda,
                 float *X, const int32_t incX);
void cblas_stbmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const int32_t K, const float *A, const int32_t lda,
                 float *X, const int32_t incX);
void cblas_stpmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const float *Ap, float *X, const int32_t incX);
void cblas_strsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const float *A, const int32_t lda, float *X,
                 const int32_t incX);
void cblas_stbsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const int32_t K, const float *A, const int32_t lda,
                 float *X, const int32_t incX);
void cblas_stpsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const float *Ap, float *X, const int32_t incX);

void cblas_dgemv(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int32_t M, const int32_t N,
                 const double alpha, const double *A, const int32_t lda,
                 const double *X, const int32_t incX, const double beta,
                 double *Y, const int32_t incY);
void cblas_dgbmv(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int32_t M, const int32_t N,
                 const int32_t KL, const int32_t KU, const double alpha,
                 const double *A, const int32_t lda, const double *X,
                 const int32_t incX, const double beta, double *Y, const int32_t incY);
void cblas_dtrmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const double *A, const int32_t lda,
                 double *X, const int32_t incX);
void cblas_dtbmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const int32_t K, const double *A, const int32_t lda,
                 double *X, const int32_t incX);
void cblas_dtpmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const double *Ap, double *X, const int32_t incX);
void cblas_dtrsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const double *A, const int32_t lda, double *X,
                 const int32_t incX);
void cblas_dtbsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const int32_t K, const double *A, const int32_t lda,
                 double *X, const int32_t incX);
void cblas_dtpsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const double *Ap, double *X, const int32_t incX);

void cblas_cgemv(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int32_t M, const int32_t N,
                 const void *alpha, const void *A, const int32_t lda,
                 const void *X, const int32_t incX, const void *beta,
                 void *Y, const int32_t incY);
void cblas_cgbmv(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int32_t M, const int32_t N,
                 const int32_t KL, const int32_t KU, const void *alpha,
                 const void *A, const int32_t lda, const void *X,
                 const int32_t incX, const void *beta, void *Y, const int32_t incY);
void cblas_ctrmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const void *A, const int32_t lda,
                 void *X, const int32_t incX);
void cblas_ctbmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const int32_t K, const void *A, const int32_t lda,
                 void *X, const int32_t incX);
void cblas_ctpmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const void *Ap, void *X, const int32_t incX);
void cblas_ctrsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const void *A, const int32_t lda, void *X,
                 const int32_t incX);
void cblas_ctbsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const int32_t K, const void *A, const int32_t lda,
                 void *X, const int32_t incX);
void cblas_ctpsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const void *Ap, void *X, const int32_t incX);

void cblas_zgemv(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int32_t M, const int32_t N,
                 const void *alpha, const void *A, const int32_t lda,
                 const void *X, const int32_t incX, const void *beta,
                 void *Y, const int32_t incY);
void cblas_zgbmv(CBLAS_LAYOUT layout,
                 CBLAS_TRANSPOSE TransA, const int32_t M, const int32_t N,
                 const int32_t KL, const int32_t KU, const void *alpha,
                 const void *A, const int32_t lda, const void *X,
                 const int32_t incX, const void *beta, void *Y, const int32_t incY);
void cblas_ztrmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const void *A, const int32_t lda,
                 void *X, const int32_t incX);
void cblas_ztbmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const int32_t K, const void *A, const int32_t lda,
                 void *X, const int32_t incX);
void cblas_ztpmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const void *Ap, void *X, const int32_t incX);
void cblas_ztrsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const void *A, const int32_t lda, void *X,
                 const int32_t incX);
void cblas_ztbsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const int32_t K, const void *A, const int32_t lda,
                 void *X, const int32_t incX);
void cblas_ztpsv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                 const int32_t N, const void *Ap, void *X, const int32_t incX);





void cblas_ssymv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const float alpha, const float *A,
                 const int32_t lda, const float *X, const int32_t incX,
                 const float beta, float *Y, const int32_t incY);
void cblas_ssbmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const int32_t K, const float alpha, const float *A,
                 const int32_t lda, const float *X, const int32_t incX,
                 const float beta, float *Y, const int32_t incY);
void cblas_sspmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const float alpha, const float *Ap,
                 const float *X, const int32_t incX,
                 const float beta, float *Y, const int32_t incY);
void cblas_sger(CBLAS_LAYOUT layout, const int32_t M, const int32_t N,
                const float alpha, const float *X, const int32_t incX,
                const float *Y, const int32_t incY, float *A, const int32_t lda);
void cblas_ssyr(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const float alpha, const float *X,
                const int32_t incX, float *A, const int32_t lda);
void cblas_sspr(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const float alpha, const float *X,
                const int32_t incX, float *Ap);
void cblas_ssyr2(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const float alpha, const float *X,
                const int32_t incX, const float *Y, const int32_t incY, float *A,
                const int32_t lda);
void cblas_sspr2(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const float alpha, const float *X,
                const int32_t incX, const float *Y, const int32_t incY, float *A);

void cblas_dsymv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const double alpha, const double *A,
                 const int32_t lda, const double *X, const int32_t incX,
                 const double beta, double *Y, const int32_t incY);
void cblas_dsbmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const int32_t K, const double alpha, const double *A,
                 const int32_t lda, const double *X, const int32_t incX,
                 const double beta, double *Y, const int32_t incY);
void cblas_dspmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const double alpha, const double *Ap,
                 const double *X, const int32_t incX,
                 const double beta, double *Y, const int32_t incY);
void cblas_dger(CBLAS_LAYOUT layout, const int32_t M, const int32_t N,
                const double alpha, const double *X, const int32_t incX,
                const double *Y, const int32_t incY, double *A, const int32_t lda);
void cblas_dsyr(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const double alpha, const double *X,
                const int32_t incX, double *A, const int32_t lda);
void cblas_dspr(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const double alpha, const double *X,
                const int32_t incX, double *Ap);
void cblas_dsyr2(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const double alpha, const double *X,
                const int32_t incX, const double *Y, const int32_t incY, double *A,
                const int32_t lda);
void cblas_dspr2(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const double alpha, const double *X,
                const int32_t incX, const double *Y, const int32_t incY, double *A);





void cblas_chemv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const void *alpha, const void *A,
                 const int32_t lda, const void *X, const int32_t incX,
                 const void *beta, void *Y, const int32_t incY);
void cblas_chbmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const int32_t K, const void *alpha, const void *A,
                 const int32_t lda, const void *X, const int32_t incX,
                 const void *beta, void *Y, const int32_t incY);
void cblas_chpmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const void *alpha, const void *Ap,
                 const void *X, const int32_t incX,
                 const void *beta, void *Y, const int32_t incY);
void cblas_cgeru(CBLAS_LAYOUT layout, const int32_t M, const int32_t N,
                 const void *alpha, const void *X, const int32_t incX,
                 const void *Y, const int32_t incY, void *A, const int32_t lda);
void cblas_cgerc(CBLAS_LAYOUT layout, const int32_t M, const int32_t N,
                 const void *alpha, const void *X, const int32_t incX,
                 const void *Y, const int32_t incY, void *A, const int32_t lda);
void cblas_cher(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const float alpha, const void *X, const int32_t incX,
                void *A, const int32_t lda);
void cblas_chpr(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const float alpha, const void *X,
                const int32_t incX, void *A);
void cblas_cher2(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, const int32_t N,
                const void *alpha, const void *X, const int32_t incX,
                const void *Y, const int32_t incY, void *A, const int32_t lda);
void cblas_chpr2(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, const int32_t N,
                const void *alpha, const void *X, const int32_t incX,
                const void *Y, const int32_t incY, void *Ap);

void cblas_zhemv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const void *alpha, const void *A,
                 const int32_t lda, const void *X, const int32_t incX,
                 const void *beta, void *Y, const int32_t incY);
void cblas_zhbmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const int32_t K, const void *alpha, const void *A,
                 const int32_t lda, const void *X, const int32_t incX,
                 const void *beta, void *Y, const int32_t incY);
void cblas_zhpmv(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 const int32_t N, const void *alpha, const void *Ap,
                 const void *X, const int32_t incX,
                 const void *beta, void *Y, const int32_t incY);
void cblas_zgeru(CBLAS_LAYOUT layout, const int32_t M, const int32_t N,
                 const void *alpha, const void *X, const int32_t incX,
                 const void *Y, const int32_t incY, void *A, const int32_t lda);
void cblas_zgerc(CBLAS_LAYOUT layout, const int32_t M, const int32_t N,
                 const void *alpha, const void *X, const int32_t incX,
                 const void *Y, const int32_t incY, void *A, const int32_t lda);
void cblas_zher(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const double alpha, const void *X, const int32_t incX,
                void *A, const int32_t lda);
void cblas_zhpr(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                const int32_t N, const double alpha, const void *X,
                const int32_t incX, void *A);
void cblas_zher2(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, const int32_t N,
                const void *alpha, const void *X, const int32_t incX,
                const void *Y, const int32_t incY, void *A, const int32_t lda);
void cblas_zhpr2(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo, const int32_t N,
                const void *alpha, const void *X, const int32_t incX,
                const void *Y, const int32_t incY, void *Ap);
# 470 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas.h"
void cblas_sgemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int32_t M, const int32_t N,
                 const int32_t K, const float alpha, const float *A,
                 const int32_t lda, const float *B, const int32_t ldb,
                 const float beta, float *C, const int32_t ldc);
void cblas_ssymm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int32_t M, const int32_t N,
                 const float alpha, const float *A, const int32_t lda,
                 const float *B, const int32_t ldb, const float beta,
                 float *C, const int32_t ldc);
void cblas_ssyrk(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                 const float alpha, const float *A, const int32_t lda,
                 const float beta, float *C, const int32_t ldc);
void cblas_ssyr2k(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                  const float alpha, const float *A, const int32_t lda,
                  const float *B, const int32_t ldb, const float beta,
                  float *C, const int32_t ldc);
void cblas_strmm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int32_t M, const int32_t N,
                 const float alpha, const float *A, const int32_t lda,
                 float *B, const int32_t ldb);
void cblas_strsm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int32_t M, const int32_t N,
                 const float alpha, const float *A, const int32_t lda,
                 float *B, const int32_t ldb);

void cblas_dgemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int32_t M, const int32_t N,
                 const int32_t K, const double alpha, const double *A,
                 const int32_t lda, const double *B, const int32_t ldb,
                 const double beta, double *C, const int32_t ldc);
void cblas_dsymm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int32_t M, const int32_t N,
                 const double alpha, const double *A, const int32_t lda,
                 const double *B, const int32_t ldb, const double beta,
                 double *C, const int32_t ldc);
void cblas_dsyrk(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                 const double alpha, const double *A, const int32_t lda,
                 const double beta, double *C, const int32_t ldc);
void cblas_dsyr2k(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                  const double alpha, const double *A, const int32_t lda,
                  const double *B, const int32_t ldb, const double beta,
                  double *C, const int32_t ldc);
void cblas_dtrmm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int32_t M, const int32_t N,
                 const double alpha, const double *A, const int32_t lda,
                 double *B, const int32_t ldb);
void cblas_dtrsm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int32_t M, const int32_t N,
                 const double alpha, const double *A, const int32_t lda,
                 double *B, const int32_t ldb);

void cblas_cgemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int32_t M, const int32_t N,
                 const int32_t K, const void *alpha, const void *A,
                 const int32_t lda, const void *B, const int32_t ldb,
                 const void *beta, void *C, const int32_t ldc);
void cblas_csymm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int32_t M, const int32_t N,
                 const void *alpha, const void *A, const int32_t lda,
                 const void *B, const int32_t ldb, const void *beta,
                 void *C, const int32_t ldc);
void cblas_csyrk(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                 const void *alpha, const void *A, const int32_t lda,
                 const void *beta, void *C, const int32_t ldc);
void cblas_csyr2k(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                  const void *alpha, const void *A, const int32_t lda,
                  const void *B, const int32_t ldb, const void *beta,
                  void *C, const int32_t ldc);
void cblas_ctrmm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int32_t M, const int32_t N,
                 const void *alpha, const void *A, const int32_t lda,
                 void *B, const int32_t ldb);
void cblas_ctrsm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int32_t M, const int32_t N,
                 const void *alpha, const void *A, const int32_t lda,
                 void *B, const int32_t ldb);

void cblas_zgemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, const int32_t M, const int32_t N,
                 const int32_t K, const void *alpha, const void *A,
                 const int32_t lda, const void *B, const int32_t ldb,
                 const void *beta, void *C, const int32_t ldc);
void cblas_zsymm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int32_t M, const int32_t N,
                 const void *alpha, const void *A, const int32_t lda,
                 const void *B, const int32_t ldb, const void *beta,
                 void *C, const int32_t ldc);
void cblas_zsyrk(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                 const void *alpha, const void *A, const int32_t lda,
                 const void *beta, void *C, const int32_t ldc);
void cblas_zsyr2k(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                  const void *alpha, const void *A, const int32_t lda,
                  const void *B, const int32_t ldb, const void *beta,
                  void *C, const int32_t ldc);
void cblas_ztrmm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int32_t M, const int32_t N,
                 const void *alpha, const void *A, const int32_t lda,
                 void *B, const int32_t ldb);
void cblas_ztrsm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, const int32_t M, const int32_t N,
                 const void *alpha, const void *A, const int32_t lda,
                 void *B, const int32_t ldb);





void cblas_chemm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int32_t M, const int32_t N,
                 const void *alpha, const void *A, const int32_t lda,
                 const void *B, const int32_t ldb, const void *beta,
                 void *C, const int32_t ldc);
void cblas_cherk(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                 const float alpha, const void *A, const int32_t lda,
                 const float beta, void *C, const int32_t ldc);
void cblas_cher2k(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                  const void *alpha, const void *A, const int32_t lda,
                  const void *B, const int32_t ldb, const float beta,
                  void *C, const int32_t ldc);

void cblas_zhemm(CBLAS_LAYOUT layout, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, const int32_t M, const int32_t N,
                 const void *alpha, const void *A, const int32_t lda,
                 const void *B, const int32_t ldb, const void *beta,
                 void *C, const int32_t ldc);
void cblas_zherk(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                 const double alpha, const void *A, const int32_t lda,
                 const double beta, void *C, const int32_t ldc);
void cblas_zher2k(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, const int32_t N, const int32_t K,
                  const void *alpha, const void *A, const int32_t lda,
                  const void *B, const int32_t ldb, const double beta,
                  void *C, const int32_t ldc);

void



cblas_xerbla(int32_t p, const char *rout, const char *form, ...);
# 11 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/src/cblas_dger.c" 2
# 1 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas_f77.h" 1
# 12 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas_f77.h"
# 1 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stdarg.h" 1 3 4
# 40 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stdarg.h" 3 4

# 40 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stdarg.h" 3 4
typedef __builtin_va_list __gnuc_va_list;
# 103 "/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-12.3.0/gcc-14.2.0-vzd2a56/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include/stdarg.h" 3 4
typedef __gnuc_va_list va_list;
# 13 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas_f77.h" 2
# 566 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas_f77.h"

# 566 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/include/cblas_f77.h"
/* xerbla_ declaration removed to avoid Tapenade pointer-analysis crash */







void srot_(const int32_t *, float *, const int32_t *, float *, const int32_t *, const float *, const float *);
void srotg_(float *,float *,float *,float *);
void srotm_(const int32_t *, float *, const int32_t *, float *, const int32_t *, const float *);
void srotmg_(float *,float *,float *,const float *, float *);
void sswap_(const int32_t *, float *, const int32_t *, float *, const int32_t *);
void scopy_(const int32_t *, const float *, const int32_t *, float *, const int32_t *);
void saxpy_(const int32_t *, const float *, const float *, const int32_t *, float *, const int32_t *);
void sdotsub_(const int32_t *, const float *, const int32_t *, const float *, const int32_t *, float *);
void sdsdotsub_(const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, float *);
void sscal_(const int32_t *, const float *, float *, const int32_t *);
void snrm2sub_(const int32_t *, const float *, const int32_t *, float *);
void sasumsub_(const int32_t *, const float *, const int32_t *, float *);
void isamaxsub_(const int32_t *, const float * , const int32_t *, int32_t *);



void drot_(const int32_t *, double *, const int32_t *, double *, const int32_t *, const double *, const double *);
void drotg_(double *,double *,double *,double *);
void drotm_(const int32_t *, double *, const int32_t *, double *, const int32_t *, const double *);
void drotmg_(double *,double *,double *,const double *, double *);
void dswap_(const int32_t *, double *, const int32_t *, double *, const int32_t *);
void dcopy_(const int32_t *, const double *, const int32_t *, double *, const int32_t *);
void daxpy_(const int32_t *, const double *, const double *, const int32_t *, double *, const int32_t *);
void dswap_(const int32_t *, double *, const int32_t *, double *, const int32_t *);
void dsdotsub_(const int32_t *, const float *, const int32_t *, const float *, const int32_t *, double *);
void ddotsub_(const int32_t *, const double *, const int32_t *, const double *, const int32_t *, double *);
void dscal_(const int32_t *, const double *, double *, const int32_t *);
void dnrm2sub_(const int32_t *, const double *, const int32_t *, double *);
void dasumsub_(const int32_t *, const double *, const int32_t *, double *);
void idamaxsub_(const int32_t *, const double * , const int32_t *, int32_t *);



void crotg_(void *, void *, float *, void *);
void csrot_(const int32_t *, void *X, const int32_t *, void *, const int32_t *, const float *, const float *);
void cswap_(const int32_t *, void *, const int32_t *, void *, const int32_t *);
void ccopy_(const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void caxpy_(const int32_t *, const void *, const void *, const int32_t *, void *, const int32_t *);
void cswap_(const int32_t *, void *, const int32_t *, void *, const int32_t *);
void cdotcsub_(const int32_t *, const void *, const int32_t *, const void *, const int32_t *, void *);
void cdotusub_(const int32_t *, const void *, const int32_t *, const void *, const int32_t *, void *);
void cscal_(const int32_t *, const void *, void *, const int32_t *);
void icamaxsub_(const int32_t *, const void *, const int32_t *, int32_t *);
void csscal_(const int32_t *, const float *, void *, const int32_t *);
void scnrm2sub_(const int32_t *, const void *, const int32_t *, float *);
void scasumsub_(const int32_t *, const void *, const int32_t *, float *);
void scabs1sub_(const void *, float *);



void zrotg_(void *, void *, double *, void *);
void zdrot_(const int32_t *, void *X, const int32_t *, void *, const int32_t *, const double *, const double *);
void zswap_(const int32_t *, void *, const int32_t *, void *, const int32_t *);
void zcopy_(const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void zaxpy_(const int32_t *, const void *, const void *, const int32_t *, void *, const int32_t *);
void zswap_(const int32_t *, void *, const int32_t *, void *, const int32_t *);
void zdotcsub_(const int32_t *, const void *, const int32_t *, const void *, const int32_t *, void *);
void zdotusub_(const int32_t *, const void *, const int32_t *, const void *, const int32_t *, void *);
void zdscal_(const int32_t *, const double *, void *, const int32_t *);
void zscal_(const int32_t *, const void *, void *, const int32_t *);
void dznrm2sub_(const int32_t *, const void *, const int32_t *, double *);
void dzasumsub_(const int32_t *, const void *, const int32_t *, double *);
void izamaxsub_(const int32_t *, const void *, const int32_t *, int32_t *);
void dcabs1sub_(const void *, double *);







void sgemv_(char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void sgbmv_(char *, const int32_t *, const int32_t *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void ssymv_(char *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void ssbmv_(char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void sspmv_(char *, const int32_t *, const float *, const float *, const float *, const int32_t *, const float *, float *, const int32_t *);
void strmv_(char *, char *, char *, const int32_t *, const float *, const int32_t *, float *, const int32_t *);
void stbmv_(char *, char *, char *, const int32_t *, const int32_t *, const float *, const int32_t *, float *, const int32_t *);
void strsv_(char *, char *, char *, const int32_t *, const float *, const int32_t *, float *, const int32_t *);
void stbsv_(char *, char *, char *, const int32_t *, const int32_t *, const float *, const int32_t *, float *, const int32_t *);
void stpmv_(char *, char *, char *, const int32_t *, const float *, float *, const int32_t *);
void stpsv_(char *, char *, char *, const int32_t *, const float *, float *, const int32_t *);
void sger_(const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, float *, const int32_t *);
void ssyr_(char *, const int32_t *, const float *, const float *, const int32_t *, float *, const int32_t *);
void sspr_(char *, const int32_t *, const float *, const float *, const int32_t *, float *);
void sspr2_(char *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, float *);
void ssyr2_(char *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, float *, const int32_t *);



void dgemv_(char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void dgbmv_(char *, const int32_t *, const int32_t *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void dsymv_(char *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void dsbmv_(char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void dspmv_(char *, const int32_t *, const double *, const double *, const double *, const int32_t *, const double *, double *, const int32_t *);
void dtrmv_(char *, char *, char *, const int32_t *, const double *, const int32_t *, double *, const int32_t *);
void dtbmv_(char *, char *, char *, const int32_t *, const int32_t *, const double *, const int32_t *, double *, const int32_t *);
void dtrsv_(char *, char *, char *, const int32_t *, const double *, const int32_t *, double *, const int32_t *);
void dtbsv_(char *, char *, char *, const int32_t *, const int32_t *, const double *, const int32_t *, double *, const int32_t *);
void dtpmv_(char *, char *, char *, const int32_t *, const double *, double *, const int32_t *);
void dtpsv_(char *, char *, char *, const int32_t *, const double *, double *, const int32_t *);
void dger_(const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, double *, const int32_t *);
void dsyr_(char *, const int32_t *, const double *, const double *, const int32_t *, double *, const int32_t *);
void dspr_(char *, const int32_t *, const double *, const double *, const int32_t *, double *);
void dspr2_(char *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, double *);
void dsyr2_(char *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, double *, const int32_t *);



void cgemv_(char *, const int32_t *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, const void *, void *, const int32_t *);
void cgbmv_(char *, const int32_t *, const int32_t *, const int32_t *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, const void *, void *, const int32_t *);
void chemv_(char *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, const void *, void *, const int32_t *);
void chbmv_(char *, const int32_t *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, const void *, void *, const int32_t *);
void chpmv_(char *, const int32_t *, const void *, const void *, const void *, const int32_t *, const void *, void *, const int32_t *);
void ctrmv_(char *, char *, char *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void ctbmv_(char *, char *, char *, const int32_t *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void ctpmv_(char *, char *, char *, const int32_t *, const void *, void *, const int32_t *);
void ctrsv_(char *, char *, char *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void ctbsv_(char *, char *, char *, const int32_t *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void ctpsv_(char *, char *, char *, const int32_t *, const void *, void *,const int32_t *);
void cgerc_(const int32_t *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void cgeru_(const int32_t *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void cher_(char *, const int32_t *, const float *, const void *, const int32_t *, void *, const int32_t *);
void cher2_(char *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void chpr_(char *, const int32_t *, const float *, const void *, const int32_t *, void *);
void chpr2_(char *, const int32_t *, const float *, const void *, const int32_t *, const void *, const int32_t *, void *);



void zgemv_(char *, const int32_t *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, const void *, void *, const int32_t *);
void zgbmv_(char *, const int32_t *, const int32_t *, const int32_t *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, const void *, void *, const int32_t *);
void zhemv_(char *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, const void *, void *, const int32_t *);
void zhbmv_(char *, const int32_t *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, const void *, void *, const int32_t *);
void zhpmv_(char *, const int32_t *, const void *, const void *, const void *, const int32_t *, const void *, void *, const int32_t *);
void ztrmv_(char *, char *, char *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void ztbmv_(char *, char *, char *, const int32_t *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void ztpmv_(char *, char *, char *, const int32_t *, const void *, void *, const int32_t *);
void ztrsv_(char *, char *, char *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void ztbsv_(char *, char *, char *, const int32_t *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void ztpsv_(char *, char *, char *, const int32_t *, const void *, void *,const int32_t *);
void zgerc_(const int32_t *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void zgeru_(const int32_t *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void zher_(char *, const int32_t *, const double *, const void *, const int32_t *, void *, const int32_t *);
void zher2_(char *, const int32_t *, const void *, const void *, const int32_t *, const void *, const int32_t *, void *, const int32_t *);
void zhpr_(char *, const int32_t *, const double *, const void *, const int32_t *, void *);
void zhpr2_(char *, const int32_t *, const double *, const void *, const int32_t *, const void *, const int32_t *, void *);







void sgemm_(char *, char *, const int32_t *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void ssymm_(char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void ssyrk_(char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, float *, const int32_t *);
void ssyr2k_(char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void strmm_(char *, char *, char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, float *, const int32_t *);
void strsm_(char *, char *, char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, float *, const int32_t *);



void dgemm_(char *, char *, const int32_t *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void dsymm_(char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void dsyrk_(char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, double *, const int32_t *);
void dsyr2k_(char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void dtrmm_(char *, char *, char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, double *, const int32_t *);
void dtrsm_(char *, char *, char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, double *, const int32_t *);



void cgemm_(char *, char *, const int32_t *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void csymm_(char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void chemm_(char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void csyrk_(char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, float *, const int32_t *);
void cherk_(char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, float *, const int32_t *);
void csyr2k_(char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void cher2k_(char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, const float *, const int32_t *, const float *, float *, const int32_t *);
void ctrmm_(char *, char *, char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, float *, const int32_t *);
void ctrsm_(char *, char *, char *, char *, const int32_t *, const int32_t *, const float *, const float *, const int32_t *, float *, const int32_t *);



void zgemm_(char *, char *, const int32_t *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void zsymm_(char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void zhemm_(char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void zsyrk_(char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, double *, const int32_t *);
void zherk_(char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, double *, const int32_t *);
void zsyr2k_(char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void zher2k_(char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, const double *, const int32_t *, const double *, double *, const int32_t *);
void ztrmm_(char *, char *, char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, double *, const int32_t *);
void ztrsm_(char *, char *, char *, char *, const int32_t *, const int32_t *, const double *, const double *, const int32_t *, double *, const int32_t *);
# 12 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/src/cblas_dger.c" 2
void cblas_dger(const CBLAS_LAYOUT layout, const int32_t M, const int32_t N,
                const double alpha, const double *X, const int32_t incX,
                const double *Y, const int32_t incY, double *A, const int32_t lda)
{

   int32_t F77_M=M, F77_N=N, F77_lda=lda, F77_incX=incX, F77_incY=incY;
# 26 "/gpfs/fs1/home/snarayan/difflapack/lapack-3.12.0-gfortran/CBLAS/src/cblas_dger.c"
   extern int CBLAS_CallFromC;
   extern int RowMajorStrg;
   RowMajorStrg = 0;

   CBLAS_CallFromC = 1;
   if (layout == CblasColMajor)
   {
      dger_(&F77_M, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A, &F77_lda)
                               ;
   }
   else if (layout == CblasRowMajor)
   {
      RowMajorStrg = 1;
      dger_(&F77_N, &F77_M ,&alpha, Y, &F77_incY, X, &F77_incX, A, &F77_lda)
                               ;

   }
   else cblas_xerbla(1, "cblas_dger", "Illegal layout setting, %d\n", layout);
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}
