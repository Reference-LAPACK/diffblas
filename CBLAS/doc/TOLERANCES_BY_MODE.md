# Step size (h) and tolerances (atol, rtol) by test mode

Generated from `run_tapenade_cblas.py` test generators. Pass criterion: `|error| <= atol + rtol * |reference|`.

## Consistency by category

| Mode | Consistent? | Exceptions |
|------|-------------|------------|
| **\_d** | **Yes** | All _d tests use the same h, atol, rtol (single: 1e-3, 2e-3, 2e-3; double: 1e-6, 1e-5, 1e-5). |
| **\_dv** | **Yes** | All _dv tests use the same values (single: 1e-3, 5e-3, 5e-3; double: 1e-6, 1e-5, 1e-5). |
| **\_b** | **No** | **nrm2** (dnrm2_b, snrm2_b): single atol=rtol=**2.0e-3**; all others: single atol=rtol=1.0e-2. Double and h are the same. |
| **\_bv** | **No** | **Scalar-result** (dasum_bv, sasum_bv, ddot_bv, sdot_bv, dnrm2_bv, snrm2_bv): double **h=1.0e-6** and single atol=rtol=**2.0e-3**; generic/gemm _bv: double h=1.0e-7, single atol=rtol=1.0e-2. |

## Full table

| Mode | Category | Precision | h (step size) | atol | rtol |
|------|----------|-----------|---------------|------|------|
| **\_d** | Forward scalar | Single (s, c) | 1.0e-3 | 2.0e-3 | 2.0e-3 |
| **\_d** | Forward scalar | Double (d, z) | 1.0e-6 | 1.0e-5 | 1.0e-5 |
| **\_dv** | Forward vector | Single (s, c) | 1.0e-3 | 5.0e-3 | 5.0e-3 |
| **\_dv** | Forward vector | Double (d, z) | 1.0e-6 | 1.0e-5 | 1.0e-5 |
| **\_b** | Reverse scalar (generic, gemm) | Single (s, c) | 1.0e-3 | 1.0e-2 | 1.0e-2 |
| **\_b** | Reverse scalar (generic, gemm) | Double (d, z) | 1.0e-7 | 1.0e-5 | 1.0e-5 |
| **\_b** | Reverse scalar (nrm2 only) | Single (s, c) | 1.0e-3 | **2.0e-3** | **2.0e-3** |
| **\_b** | Reverse scalar (nrm2 only) | Double (d, z) | 1.0e-7 | 1.0e-5 | 1.0e-5 |
| **\_bv** | Reverse vector (generic VJP, gemm) | Single (s, c) | 1.0e-3 | 1.0e-2 | 1.0e-2 |
| **\_bv** | Reverse vector (generic VJP, gemm) | Double (d, z) | 1.0e-7 | 1.0e-5 | 1.0e-5 |
| **\_bv** | Reverse vector (scalar-result: dasum, ddot, nrm2, etc.) | Single (s, c) | 1.0e-3 | **2.0e-3** | **2.0e-3** |
| **\_bv** | Reverse vector (scalar-result) | Double (d, z) | **1.0e-6** | 1.0e-5 | 1.0e-5 |

## Notes

- **\_d**: Matches Fortran BLAS forward tests (e.g. `test_sgemm.f90` / `test_dgemm.f90`).
- **\_dv**: Same h as _d; atol/rtol slightly looser for single precision (5.0e-3) for multi-direction FD.
- **\_b** / **\_bv** (generic): VJP check; smaller h (1.0e-7 for double) for central-difference stability; looser single-precision tolerances (1.0e-2).
- **\_bv** (nrm2-style): Used for scalar-result routines (e.g. snrm2_bv, dnrm2_bv); h and atol/rtol aligned with nrm2 _b/_dv tests.
- **nrm2 _b** (reverse scalar): Same as generic _b (h=1.0e-7 double / 1.0e-3 float; atol=rtol=1.0e-5 double / 1.0e-2 float). **nrm2 _d** uses h=1.0e-7 double, atol=rtol=1.0e-5; single uses h=1.0e-3, atol=rtol=2.0e-3.
