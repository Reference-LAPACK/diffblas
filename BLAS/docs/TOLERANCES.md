# Differentiation test tolerances

Tolerances and step sizes for finite-difference derivative checks in the BLAS differentiation test generator.

---

## Defaults

### rtol/atol by precision family

| Family | Meaning | rtol | atol |
|--------|---------|------|------|
| S      | `S*` (single real) | 2.0e-3 | 2.0e-3 |
| C      | `C*` (single complex) | 1.0e-3 | 1.0e-3 |
| D      | `D*` (double real) | 1.0e-5 | 1.0e-5 |
| Z      | `Z*` (double complex) | 1.0e-5 | 1.0e-5 |

### step size h by precision family

| Family | h |
|--------|---|
| S, C   | 1.0e-3 |
| D, Z   | 1.0e-7 |

---

## Overrides

### Mixed-precision D* (single-precision first differentiable input)

Applies when the routine behaves like “double output, but first differentiable input is single precision” (e.g. `DSDOT` with **SX** first; the generator also treats **SY** and **SB** as single-precision inputs for `D*`).

- **Scalar forward**: override **h = 1.0e-3** (rtol/atol remain `D*` base = 1.0e-5)
- **Scalar reverse / vector forward / vector reverse**: override **h = 1.0e-3**, **rtol = atol = 2.0e-3**

### Relaxed C* tolerance in vector reverse

Only for **single-precision complex** (`C*`) **vector reverse** tests:

| Routine family (examples) | rtol/atol |
|---------------------------|-----------|
| DOT (e.g. `CDOTC`)        | 2.5e-2 |
| BLAS3 (e.g. `CGEMM`, `CSYMM`, `CHEMM`) | 1.0e-2 |
| BLAS2 banded MV (e.g. `CGBMV`, `CTBMV`, `CHBMV`) | 1.0e-2 |

All other `C*` modes use the base tolerance (1.0e-3). `Z*` does not use relaxed tolerances.
