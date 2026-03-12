# Differentiation test tolerances

Tolerances and step sizes used for finite-difference checks in BLAS differentiation tests (scalar/vector, forward/reverse). All modes use the same precision-based scheme unless a mixed-precision override applies.

---

## Base tolerances by precision

| Family | Description           | rtol    | atol    |
|--------|-----------------------|---------|---------|
| S      | single real (`S*`)    | 2.0e-3  | 2.0e-3  |
| C      | single complex (`C*`) | 1.0e-3  | 1.0e-3  |
| D      | double real (`D*`)    | 1.0e-5  | 1.0e-5  |
| Z      | double complex (`Z*`) | 1.0e-5  | 1.0e-5  |

These values are used in:

- Scalar forward
- Scalar reverse
- Vector forward
- Vector reverse

---

## Step size (h)

For non–mixed-precision functions:

| Precision   | h        |
|------------|----------|
| S*, C*     | 1.0e-3   |
| D*, Z*     | 1.0e-7   |

(≈ 10·√ε for double precision.)

---

## Mixed-precision override

For routines whose **output is double precision** but whose **first differentiable input** is **single precision** (e.g. `DSDOT`), the generator uses single-precision–style settings so the finite-difference check matches the conditioning of the inputs:

- **h** = 1.0e-3  
- **rtol** = 2.0e-3  
- **atol** = 2.0e-3  

This override is applied in:

- Scalar reverse
- Vector forward
- Vector reverse  

Detection: `precision_type == real(8)` and the first entry in the `inputs` list has `get_param_precision(first_input, func_name, param_types) == "real(4)"`. In the generator, `get_param_precision` returns `real(4)` for **D\*** functions when the parameter is one of **SX**, **SY**, **SB**.

---

## Mixed-precision tests (list)

A test is treated as mixed-precision if it is for a **D\*** (or **Z\***) routine and the **first differentiable input** is single precision. The generator explicitly treats **SX**, **SY**, and **SB** as single precision for **D\*** routines.

**Routines that use the mixed-precision override** (when present in the suite and documented with that input order):

| Routine | First input(s) | Modes using override        |
|---------|----------------|-----------------------------|
| **DSDOT** | SX (then SY) | Scalar reverse, vector forward, vector reverse |

**Note:** Any other **D\*** routine whose first `\param[in]` is **SX**, **SY**, or **SB** will also get the override. There is no **Z\*** branch for single-precision inputs in `get_param_precision`, so currently only **D\*** routines can be mixed-precision in this sense. If you add a **D\*** (or in future **Z\***) routine with a single-precision first input, it will automatically receive the same h and tolerances as above.

---

## Summary table (all modes)

| Mode             | S* / C* (h) | D* / Z* (h) | Mixed-precision (h, rtol, atol)      |
|------------------|-------------|-------------|---------------------------------------|
| Scalar forward   | 1e-3 / 2e-3 or 1e-3 | 1e-7 / 1e-5 | h = 1e-3 only (rtol/atol stay 1e-5)   |
| Scalar reverse   | 1e-3 / 2e-3 or 1e-3 | 1e-7 / 1e-5 | 1e-3, 2e-3, 2e-3                      |
| Vector forward   | 1e-3 / 2e-3 or 1e-3 | 1e-7 / 1e-5 | 1e-3, 2e-3, 2e-3                      |
| Vector reverse   | 1e-3 / 2e-3 or 1e-3 | 1e-7 / 1e-5 | 1e-3, 2e-3, 2e-3                      |

(Base tolerances for S/C/D/Z are as in the first table; mixed-precision replaces h and rtol/atol only where indicated. In scalar forward, mixed-precision only changes the step size h to 1e-3; rtol/atol remain 1e-5.)
