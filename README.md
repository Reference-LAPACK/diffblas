# `diffblas`

`diffblas`is a library that  provides (algorithmically) differentiated BLAS routines from their reference implementation in [lapack](https://github.com/Reference-LAPACK/lapack) using the automatic differentiation tool [Tapenade](https://gitlab.inria.fr/tapenade/tapenade) in four modes: forward (`_d`), vector forward (`_dv`), reverse (`_b`), and vector reverse (`_bv`).
The compiled `libdiffblas` can be linked into applications that need derivatives of BLAS operations for optimization, sensitivity analysis etc.

This work was inspired in part by a need to differentiate a Fortran code [HFBTHO](https://www.sciencedirect.com/science/article/abs/pii/S0010465525004564), that uses LAPACK and BLAS routines, and to use the differentiated application for gradient-based optimization. 

## Using the pre-compiled library from your application

Use the library the same way the tests do: 
1. Call the differentiated routines from your application.
2. Include the right header/module.
3. Link your application with `libdiffblas` and the BLAS library (`refblas`).


**Calling convention (same as in the tests):**

1. **Forward (`_d`):** Declare the differentiated routine as `external` (or use a module). Call it with the same arguments as the original BLAS routine, plus one extra “derivative” argument per floating-point input/output (e.g. `a`, `a_d`). Example for DGEMM:
   ```fortran
   external :: dgemm_d
   ! ... set transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc
   ! ... set alpha_d, a_d, b_d, beta_d, c_d (derivative inputs; c_d is output)
   call dgemm_d(transa, transb, m, n, k, alpha, alpha_d, a, a_d, lda, b, b_d, ldb, beta, beta_d, c, c_d, ldc)
   ```
2. **Reverse (`_b`):** Call the reverse routine with the same shape as the forward call; the “b” arguments carry adjoints (see generated test files such as `BLAS/test/test_dgemm_reverse.f90`).
3. **Vector forward (`_dv`) / vector reverse (`_bv`):** Same idea with an extra dimension for multiple directions; see `BLAS/test/test_*_vector_forward.f90` and `test_*_vector_reverse.f90`. 

***Handing dimensions :***
BLAS uses assumed size arrays as seen in `dgemm.f`. 
```fortran
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
```

In order to enable the use of the differentiated versions of these routines in a portable way, we require you to call addtional routines before calling the differentiated rouint. For **reverse** (`_b`) and **vector reverse** (`_bv`) — and in some cases **vector forward** (`_dv`) the differentiated code expects certain "ISIZE" values to be set before the call. These describe the dimensions of the assumed size arrays and matrices. 

For example (as in `BLAS/test/test_dgemm_reverse.f90`): for `dgemm_b`, set the second dimension of the two matrix arguments with the appropriate dimensions:

```fortran
use DIFFSIZES
! ...
call set_ISIZE2OFA(lda)   ! leading dimension of A
call set_ISIZE2OFB(ldb)   ! leading dimension of B
call dgemm_b(transa, transb, m, n, k, alpha, alphab, a, ab, lda, b, bb, ldb, beta, betab, c, cb, ldc)
call set_ISIZE2OFA(-1)
call set_ISIZE2OFB(-1)
```

Matrix routines often use `set_ISIZE2OFA` / `set_ISIZE2OFB`; level-1 with vectors may use `set_ISIZE1OFX`, `set_ISIZE1OFY`, or `set_ISIZE1OFZx` / `set_ISIZE1OFZy`). The  tests in `BLAS/test/` show the exact `set_ISIZE*` calls for each routine. If a setter is not called, the differentiated will stop with an error. You can pptionally reset the dimensions by calling set_ISIZE*(-1)` so that the next use must set the dimension again. This is recommended if you reuse the same differentiated routine in multiple contexts.

***Handing directions in vector mode:***
We have currently hardcoded the number of directions to 4. In case you need a different number, you should download the source of this repository and edit both `BLAS/include/DIFFSIZES.inc` and `BLAS/include/DIFFSIZES.F90`. Then build the library to generate the appropriate library.
```fortran
INTEGER, PARAMETER :: nbdirsmax = 4
```
```fortran
      integer nbdirsmax
      parameter (nbdirsmax=4)
```fortran


**Minimal example (forward mode, one routine):** Compile your main program and link with the built library and BLAS:

```bash
gfortran -O2 -I/path/to/diff-lapack/BLAS/include -c my_main.f90 -o my_main.o
gfortran -o my_main my_main.o /path/to/diff-lapack/builddir/libdiffblas.a -L$LAPACKDIR -lrefblas
```

(If you built with Make, use `BLAS/build/libdiffblas_d.a` for forward mode only.) The tests in `BLAS/test/` (e.g. `test_dgemm.f90`, `test_dgemm_reverse.f90`) are full examples of how to declare and call the differentiated routines and how to link; use them as templates for your application.

## Building the library (and tests)

You need **pre-generated** sources (from step 1). The build compiles them and links with a BLAS library and the Tapenade adStack runtime.

### Build with Meson (library and tests)

**Dependencies:**

- **Fortran compiler** (e.g. gfortran, ifort, ifx) and **C compiler** (e.g. gcc).
- **LAPACK installation** — a built Reference LAPACK (or compatible) providing BLAS (e.g. `librefblas.a` or `libblas.a`). Set **`LAPACKDIR`** (or equivalent) so Meson can find it (see below).
- **Tapenade adStack** — the repo already contains `TAPENADE/adStack.c` and `TAPENADE/include/`; Meson compiles and links these automatically. No separate Tapenade install is required for the build.

**Configure and build from the project root:**

```bash
# If BLAS is in a custom location (e.g. LAPACK build dir), pass library search path
meson setup builddir -Dlibblas=refblas -Dlibblas_path=/path/to/lapack/build

# Or rely on system / environment (e.g. LIBRARY_PATH, or refblas in default path)
meson setup builddir -Dlibblas=refblas

meson compile -C builddir
```

This produces `builddir/libdiffblas.a` (or shared library if configured with `-Ddefault_library=shared`). The Meson build uses `BLAS/meson.build`, compiles everything in `BLAS/src/` plus `BLAS/include/DIFFSIZES.f90`, `BLAS/src/DIFFSIZES_access.f`, and `TAPENADE/adStack.c`; it does **not** build the test executables (tests are built by the BLAS Makefile in 2b if you use that).

To also run the tests, use the BLAS Makefile (2b) or run the test programs built there. To install the library:

```bash
meson install -C builddir --prefix /your/install
```

### Build with Make (library and tests)

**Dependencies:**

- **Fortran compiler** (e.g. gfortran) and **C compiler** (e.g. gcc).
- **LAPACK installed** — build and install [Reference LAPACK](https://github.com/Reference-LAPACK/lapack) so that you have `librefblas` (or set `BLAS_LIB` in the Makefile to point at your BLAS).
- **`LAPACKDIR`** — set to the directory where LAPACK was built (or where `librefblas.a` lives), so the linker can find `-lrefblas`:
  ```bash
  export LAPACKDIR=/path/to/lapack/build
  ```
- **Tapenade adStack** — the Makefile looks for `adStack.c` in this order: `BLAS/src/adStack.c`, then `../TAPENADE/adStack.c`, then `$TAPENADEDIR/ADFirstAidKit/adStack.c`. The repo’s `TAPENADE/` copy is used by default; no Tapenade install needed if you do not override.

**Build from the BLAS directory:**

```bash
cd BLAS
export LAPACKDIR=/path/to/your/lapack/build   # or wherever librefblas is
make
```

This builds per-mode static libraries (`build/libdiffblas_d.a`, `libdiffblas_b.a`, `libdiffblas_dv.a`, `libdiffblas_bv.a`) and test executables in `build/`. Run tests:

```bash
./run_tests.sh
# or run individual test executables, e.g. build/test_dgemm, build/test_dgemm_reverse
```

## Creating the library sources with run_tapenade_blas.py

This step **generates** the differentiated Fortran sources and test programs (e.g. under `BLAS/src/`, `BLAS/test/`, or a custom output directory like `out/`). 

**Dependencies:**

- **LAPACK source tree** — Reference LAPACK (or at least its BLAS SRC directory) so that `--input-dir` points to the BLAS source files (e.g. `lapack-3.x/BLAS/SRC/`).
- **Tapenade** — installed and on your `PATH` (or pass `--tapenade-bin=/path/to/tapenade`).
- **Python 3** — to run `run_tapenade_blas.py`.

**Example: generate from a Reference LAPACK tarball**

```bash
# Download and unpack Reference LAPACK (or set LAPACK_SRC to your tree)
wget -q -O tmp.html https://raw.githubusercontent.com/Reference-LAPACK/lapack/refs/heads/master/README.md
VERSION=$(grep VERSION tmp.html | tail -1 | awk '{print $3}')
rm tmp.html
wget -P tmp/ https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v${VERSION}.tar.gz
tar xzf tmp/v${VERSION}.tar.gz -C tmp/

# Generate differentiated BLAS into BLAS/ (flat layout: src/, test/, include/)
python run_tapenade_blas.py --input-dir tmp/lapack-${VERSION}/BLAS/SRC/ --out-dir BLAS --flat
```

**Options:** Use `--file dgemm.f sgemm.f` to restrict to specific routines; `--mode d b` to generate only certain modes; see `python run_tapenade_blas.py --help`.
  
## Future work
There are several routines within BLAS that are not differentiated with verification, which we plan to add. We will add differentiated versions of the CBLAS routines. 
We are beginning to write a differentiated version of the LAPACK library.

## Acknowledgement
This work was supported in part by the Applied Mathematics activity within the U.S. Department of Energy, Office of Science, Office
of Advanced Scientific Computing Research Applied Mathematics, and Office of Nuclear Physics SciDAC program under Contract No. DE-AC02-06CH11357. This work was supported in part by NSF CSSI grant 2104068. 
