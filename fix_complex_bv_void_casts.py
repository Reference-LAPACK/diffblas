#!/usr/bin/env python3
"""
Post-process Tapenade-generated complex _bv.c files to fix void* dereference errors.

Tapenade emits void* for the C BLAS complex API; the generated code then incorrectly
dereferences or subscripts those pointers (e.g. *alphab[nd] = 0.0), which is invalid in C.

This script replaces those patterns with proper casts:
- Single precision complex (cblas_c*_bv.c): use (float complex *) e.g. ((float complex *)alphab)[nd]
- Double precision complex (cblas_z*_bv.c): use (double complex *) e.g. ((double complex *)alphab)[nd]

Usage:
  python fix_complex_bv_void_casts.py [src_dir]
  Default src_dir: out_cblas/src (relative to script dir) or current directory.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

# Adjoint parameter names that appear as void* in Tapenade output and get dereferenced
ADJOINT_PARAMS = frozenset({"alphab", "betab", "Ab", "Bb", "Cb", "Xb", "Yb", "dotcb", "dotub", "Apb"})

# Matrix adjoint params that are declared as (*Param)[NBDirsMax] in real _bv.c (e.g. dsyrk, dsyr2k).
# Param[nd] = 0.0 is invalid (array type); replace with loop zeroing Param[nd][_ii].
# Only Ab is known to have this declaration in dsyrk/dsyr2k; Bb in dsyr2k is double* so leave it.
ARRAY_TYPE_ADJOINT_PARAMS = ("Ab",)


def fix_complex_bv_file(path: Path, is_single: bool) -> bool:
    """Apply void* cast fixes to one _bv.c file. Returns True if file was modified."""
    cast_type = "float complex" if is_single else "double complex"
    text = path.read_text()
    original = text

    # Build regex pattern for *param[nd] and (*param)[nd] for each known adjoint name.
    # Match *alphab[nd] or (*alphab)[nd] and replace with ((float complex *)alphab)[nd]
    for param in ADJOINT_PARAMS:
        # Pattern 1: *param[nd] (invalid: dereferencing void)
        # Replace with ((cast_type *)param)[nd]
        text = re.sub(
            rf"\*{re.escape(param)}\[nd\]",
            f"(({cast_type} *){param})[nd]",
            text,
        )
        # Pattern 2: (*param)[nd] (invalid: subscripting void*)
        text = re.sub(
            rf"\(\*{re.escape(param)}\)\[nd\]",
            f"(({cast_type} *){param})[nd]",
            text,
        )
        # Pattern 3: *((cast_type *)param)[nd] (redundant cast+deref; normalize to ((cast_type *)param)[nd])
        text = re.sub(
            rf"\*\(\s*\(\s*{re.escape(cast_type)}\s*\*\s*\)\s*{re.escape(param)}\s*\)\[nd\]",
            f"(({cast_type} *){param})[nd]",
            text,
        )
        # Pattern 4: param[nd] = 0.0 (bare subscript on void* - no leading * or parens; e.g. Ab[nd] = 0.0, Bb[nd] = 0.0)
        text = re.sub(
            rf"\b{re.escape(param)}\[nd\]\s*=\s*0\.0",
            f"(({cast_type} *){param})[nd] = 0.0",
            text,
        )

    # Remove commented-out incorrect Tapenade lines (e.g. //*alphab[nd] = 0.0;)
    # and duplicate comment lines that mirror the fix (e.g. //((float complex *)alphab)[nd] = 0.0;)
    for param in ADJOINT_PARAMS:
        text = re.sub(
            rf"\n\s*//\s*\*{re.escape(param)}\[nd\][^\n]*",
            "",
            text,
        )
        text = re.sub(
            rf"\n\s*//\s*\(\s*\(\s*{re.escape(cast_type)}\s*\*\s*\)\s*{re.escape(param)}\s*\)\[nd\][^\n]*",
            "",
            text,
        )

    # --- c*_bv (single precision) only: fix Tapenade's void* and local pointer-array misuse ---
    if is_single:
        # Local declarations: Tapenade uses float *xb; float *yb; (and sometimes yyb) but then xb[nd], yb[nd] as pointers -> declare as array of pointers
        text = re.sub(
            r"\bfloat\s+\*xb\s*;",
            "float *xb[NBDirsMax];",
            text,
        )
        text = re.sub(
            r"\bfloat\s+\*yb\s*;",
            "float *yb[NBDirsMax];",
            text,
        )
        text = re.sub(
            r"\bfloat\s+\*yyb\s*;",
            "float *yyb[NBDirsMax];",
            text,
        )
        # Xb[nd], Yb[nd] are invalid (void* subscript). API passes void* Xb; treat as (float *)Xb + nd for direction nd base.
        text = re.sub(
            r"xb\[nd\]\s*=\s*\(float\s*\*\)Xb\[nd\]",
            "xb[nd] = (float *)((float *)Xb + nd)",
            text,
        )
        text = re.sub(
            r"yb\[nd\]\s*=\s*\(float\s*\*\)Yb\[nd\]",
            "yb[nd] = (float *)((float *)Yb + nd)",
            text,
        )
        text = re.sub(
            r"yyb\[nd\]\s*=\s*\(float\s*\*\)Yb\[nd\]",
            "yyb[nd] = (float *)((float *)Yb + nd)",
            text,
        )
        # float *xxb = (const float *)Xb[nd] -> 2D view so xxb[i][nd] works
        text = re.sub(
            r"float\s+\*xxb\s*=\s*\(const\s+float\s*\*\)Xb\[nd\]",
            "float (*xxb)[NBDirsMax] = (float (*)[NBDirsMax])((void *)Xb)",
            text,
        )
        # xb[nd] / yb[nd] = (float (*)[NBDirsMax])malloc(...) -> (float *)malloc(n*sizeof(float[NBDirsMax]))
        for arr in ("xb", "yb"):
            text = re.sub(
                rf"{re.escape(arr)}\[nd\]\s*=\s*\(float\s*\(\*\)\[NBDirsMax\]\)malloc\s*\(\s*n\s*\*\s*sizeof\s*\(\s*float\s*\[\s*",
                f"{arr}[nd] = (float *)malloc(n*sizeof(float [",
                text,
            )
        text = re.sub(
            r"(\s*NBDirsMax\s*\]\s*\)\s*)\*\s*NBDirsMax\s*\)",
            r"\1)",
            text,
        )
        # Zeroing loop: xb[ii1][nd] / yb[ii1][nd] = 0.0 with float *xb[] -> xb[nd][ii1] = 0.0
        text = re.sub(
            r"\bxb\[ii1\]\[nd\]\s*=\s*0\.0",
            "xb[nd][ii1] = 0.0",
            text,
        )
        text = re.sub(
            r"\byb\[ii1\]\[nd\]\s*=\s*0\.0",
            "yb[nd][ii1] = 0.0",
            text,
        )
        # Parenthesized expr with nested parens e.g. (xb+(n-2))[]: fix ")[]" -> ")[nd]"
        text = re.sub(r"\)\s*\[\]", ")[nd]", text)
        # xxb[nd] = (xxb+i)[nd] is array assignment -> memcpy
        text = re.sub(
            r"xxb\[nd\]\s*=\s*\(xxb\+i\)\[nd\]",
            "memcpy(xxb[nd], (xxb+i)[nd], NBDirsMax*sizeof(float))",
            text,
        )
        # Pointer arithmetic forms xxb[nd] = (xxb+...) [nd]
        text = re.sub(
            r"xxb\[nd\]\s*=\s*\((xxb\+[^)]+)\)\[nd\]",
            r"memcpy(xxb[nd], (\1)[nd], NBDirsMax*sizeof(float))",
            text,
        )
        # xb[nd] = (xb+...)[nd] and similar - pointer assignment, RHS is float*; (xb+(n-2))[nd] etc. need cast for clarity
        text = re.sub(
            r"xb\[nd\]\s*=\s*\((xb\+[^)]+)\)\[nd\]",
            r"xb[nd] = (\1)[nd]",
            text,
        )
        text = re.sub(
            r"xb\[nd\]\s*=\s*\((xb\+tincx)\)\[nd\]",
            r"xb[nd] = (xb+tincx)[nd]",
            text,
        )
        # (*yb)[nd] and (*xb)[nd] when var is restored float* -> ((float *)yb)[nd]
        text = re.sub(
            r"\(\*yb\)\[nd\]",
            "((float *)yb)[nd]",
            text,
        )
        text = re.sub(
            r"\(\*xb\)\[nd\]",
            "((float *)xb)[nd]",
            text,
        )
        text = re.sub(
            r"\(\*xxb\)\[nd\]",
            "((float *)xxb)[nd]",
            text,
        )
        text = re.sub(
            r"\(\*yyb\)\[nd\]",
            "((float *)yyb)[nd]",
            text,
        )
        # In pop block, xb[1][nd] or yyb[1][nd] when var was restored to float* -> ((float *)var)[nd]
        text = re.sub(
            r"\bxb\[1\]\[nd\]",
            "((float *)xb)[nd]",
            text,
        )
        text = re.sub(
            r"\byyb\[1\]\[nd\]",
            "((float *)yyb)[nd]",
            text,
        )
        # Invalid LHS and RHS: betab/alphab are void* to 2*NBDirsMax floats (real,imag per direction).
        # Use (float *) + linear index so LHS is clearly scalar: [i][nd] -> ((float *)(var))[(i)*NBDirsMax+(nd)]
        def _scalar_2d(var: str, i: str, j: str) -> str:
            return f"((float *)({var}))[({i})*NBDirsMax+({j})]"
        text = re.sub(
            r"\(const\s+float\s*\*\)betab\[1\]\[nd\]",
            _scalar_2d("betab", "1", "nd"),
            text,
        )
        text = re.sub(
            r"\(\*\s*\(\s*const\s+float\s*\*\s*\)betab\)\[nd\]",
            _scalar_2d("betab", "0", "nd"),
            text,
        )
        text = re.sub(
            r"\(const\s+float\s*\*\)alphab\[1\]\[nd\]",
            _scalar_2d("alphab", "1", "nd"),
            text,
        )
        text = re.sub(
            r"\(\*\s*\(\s*const\s+float\s*\*\s*\)alphab\)\[nd\]",
            _scalar_2d("alphab", "0", "nd"),
            text,
        )
        # Normalize any ((float complex *)betab)[1][nd] / [nd] from earlier pass to the 2D layout (scalar form)
        text = re.sub(
            r"\(\(\s*float\s+complex\s+\*\s*\)betab\)\[1\]\[nd\]",
            _scalar_2d("betab", "1", "nd"),
            text,
        )
        text = re.sub(
            r"\(\(\s*float\s+complex\s+\*\s*\)betab\)\[nd\]",
            _scalar_2d("betab", "0", "nd"),
            text,
        )
        text = re.sub(
            r"\(\(\s*float\s+complex\s+\*\s*\)alphab\)\[1\]\[nd\]",
            _scalar_2d("alphab", "1", "nd"),
            text,
        )
        text = re.sub(
            r"\(\(\s*float\s+complex\s+\*\s*\)alphab\)\[nd\]",
            _scalar_2d("alphab", "0", "nd"),
            text,
        )
        # Scalar adjoint zeroing: ((float complex *)alphab)[nd] = 0.0 -> zero both real and imag (use block so LHS is scalar)
        text = re.sub(
            r"\(\(\s*float\s+complex\s+\*\s*\)alphab\)\[nd\]\s*=\s*0\.0",
            "do { " + _scalar_2d("alphab", "0", "nd") + " = 0.0; " + _scalar_2d("alphab", "1", "nd") + " = 0.0; } while(0)",
            text,
        )
        text = re.sub(
            r"\(\(\s*float\s+complex\s+\*\s*\)betab\)\[nd\]\s*=\s*0\.0",
            "do { " + _scalar_2d("betab", "0", "nd") + " = 0.0; " + _scalar_2d("betab", "1", "nd") + " = 0.0; } while(0)",
            text,
        )
        # Fix any remaining (...)[0][nd] = 0.0 that came from comma form (convert to scalar form)
        text = re.sub(
            r"\(\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*alphab\s*\)\)\[0\]\[nd\]\s*=\s*0\.0\s*,\s*\(\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*alphab\s*\)\)\[1\]\[nd\]\s*=\s*0\.0",
            "do { (((float (*)[2][NBDirsMax])(alphab))[0])[nd] = 0.0; (((float (*)[2][NBDirsMax])(alphab))[1])[nd] = 0.0; } while(0)",
            text,
        )
        text = re.sub(
            r"\(\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*betab\s*\)\)\[0\]\[nd\]\s*=\s*0\.0\s*,\s*\(\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*betab\s*\)\)\[1\]\[nd\]\s*=\s*0\.0",
            "do { (((float (*)[2][NBDirsMax])(betab))[0])[nd] = 0.0; (((float (*)[2][NBDirsMax])(betab))[1])[nd] = 0.0; } while(0)",
            text,
        )
        # Normalize to linear indexing so LHS is scalar: (...)[0][nd] or ((...)[0])[nd] -> ((float *)(var))[0*NBDirsMax+nd]
        text = re.sub(
            r"\(\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*betab\s*\)\)\[0\]\)\[nd\]",
            "((float *)(betab))[0*NBDirsMax+nd]",
            text,
        )
        text = re.sub(
            r"\(\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*betab\s*\)\)\[1\]\)\[nd\]",
            "((float *)(betab))[1*NBDirsMax+nd]",
            text,
        )
        text = re.sub(
            r"\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*betab\s*\)\)\[0\]\[nd\]",
            "((float *)(betab))[0*NBDirsMax+nd]",
            text,
        )
        text = re.sub(
            r"\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*betab\s*\)\)\[1\]\[nd\]",
            "((float *)(betab))[1*NBDirsMax+nd]",
            text,
        )
        text = re.sub(
            r"\(\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*alphab\s*\)\)\[0\]\)\[nd\]",
            "((float *)(alphab))[0*NBDirsMax+nd]",
            text,
        )
        text = re.sub(
            r"\(\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*alphab\s*\)\)\[1\]\)\[nd\]",
            "((float *)(alphab))[1*NBDirsMax+nd]",
            text,
        )
        text = re.sub(
            r"\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*alphab\s*\)\)\[0\]\[nd\]",
            "((float *)(alphab))[0*NBDirsMax+nd]",
            text,
        )
        text = re.sub(
            r"\(\(float\s*\(\*\)\[2\]\[NBDirsMax\]\)\s*\(\s*alphab\s*\)\)\[1\]\[nd\]",
            "((float *)(alphab))[1*NBDirsMax+nd]",
            text,
        )
        # yb[nd]++ and --yb[nd]: yb is float*[] so yb[nd] is float*; incrementing pointer is valid
        # (no change needed)
        # Ensure stdlib.h for memcpy if not present
        if "memcpy" in text and "#include <stdlib.h>" in text and "#include <string.h>" not in text:
            text = text.replace("#include <stdlib.h>", "#include <stdlib.h>\n#include <string.h>", 1)

    # --- z*_bv (double precision complex): same fixes with double, and Tapenade's empty [] -> [nd] ---
    else:
        # Assignment fixes first (before [] -> [nd]): xb[] = (double *)Xb[] etc.
        text = re.sub(
            r"xb\[\]\s*=\s*\(double\s*\*\)Xb\[\]",
            "xb[nd] = (double *)((double *)Xb + nd)",
            text,
        )
        text = re.sub(
            r"yb\[\]\s*=\s*\(double\s*\*\)Yb\[\]",
            "yb[nd] = (double *)((double *)Yb + nd)",
            text,
        )
        text = re.sub(
            r"yyb\[\]\s*=\s*\(double\s*\*\)Yb\[\]",
            "yyb[nd] = (double *)((double *)Yb + nd)",
            text,
        )
        # Tapenade sometimes emits xb[], yb[] (empty brackets); normalize to xb[nd], yb[nd]
        for var in ("xb", "yb", "yyb", "xxb"):
            text = re.sub(rf"\b{re.escape(var)}\[\]", f"{var}[nd]", text)
        # RHS expressions with []: (xb+...)[], (yb+...)[], txb[], tyb[], (xxb+i)[], etc.
        # Parenthesized expr with nested parens e.g. (xb+(n-2))[]: fix ")[]" -> ")[nd]" first
        text = re.sub(r"\)\s*\[\]", ")[nd]", text)
        text = re.sub(r"\((xb\+[^)]+)\)\[\]", r"(\1)[nd]", text)
        text = re.sub(r"\((yb\+[^)]+)\)\[\]", r"(\1)[nd]", text)
        text = re.sub(r"\((yb-n)\)\[\]", r"(yb-n)[nd]", text)
        text = re.sub(r"\((xb-n)\)\[\]", r"(xb-n)[nd]", text)
        text = re.sub(r"\((xxb\+[^)]+)\)\[\]", r"(\1)[nd]", text)
        text = re.sub(r"\((yyb\+[^)]+)\)\[\]", r"(\1)[nd]", text)
        text = re.sub(r"\btxb\[\]", "txb[nd]", text)
        text = re.sub(r"\btyb\[\]", "tyb[nd]", text)
        # Local declarations: double *xb; -> double *xb[NBDirsMax]; (and yb, yyb, xxb, txb, tyb)
        text = re.sub(r"\bdouble\s+\*xb\s*;", "double *xb[NBDirsMax];", text)
        text = re.sub(r"\bdouble\s+\*yb\s*;", "double *yb[NBDirsMax];", text)
        text = re.sub(r"\bdouble\s+\*yyb\s*;", "double *yyb[NBDirsMax];", text)
        text = re.sub(r"\bdouble\s+\*txb\s*;", "double *txb[NBDirsMax];", text)
        text = re.sub(r"\bdouble\s+\*tyb\s*;", "double *tyb[NBDirsMax];", text)
        text = re.sub(
            r"\bdouble\s+\*xxb\s*;",
            "double (*xxb)[NBDirsMax] = (double (*)[NBDirsMax])((void *)Xb);",
            text,
        )
        # xb[nd] = (double *)Xb[nd] when [nd] already present
        text = re.sub(
            r"xb\[nd\]\s*=\s*\(double\s*\*\)Xb\[nd\]",
            "xb[nd] = (double *)((double *)Xb + nd)",
            text,
        )
        text = re.sub(
            r"yb\[nd\]\s*=\s*\(double\s*\*\)Yb\[nd\]",
            "yb[nd] = (double *)((double *)Yb + nd)",
            text,
        )
        text = re.sub(
            r"yyb\[nd\]\s*=\s*\(double\s*\*\)Yb\[nd\]",
            "yyb[nd] = (double *)((double *)Yb + nd)",
            text,
        )
        # double *xxb = (const double *)Xb[nd] -> 2D view
        text = re.sub(
            r"double\s+\*xxb\s*=\s*\(const\s+double\s*\*\)Xb\[nd\]",
            "double (*xxb)[NBDirsMax] = (double (*)[NBDirsMax])((void *)Xb)",
            text,
        )
        # xb[nd] / yb[nd] = (double (*)[NBDirsMax])malloc(...)
        for arr in ("xb", "yb"):
            text = re.sub(
                rf"{re.escape(arr)}\[nd\]\s*=\s*\(double\s*\(\*\)\[NBDirsMax\]\)malloc\s*\(\s*n\s*\*\s*sizeof\s*\(\s*double\s*\[\s*",
                f"{arr}[nd] = (double *)malloc(n*sizeof(double [",
                text,
            )
        text = re.sub(
            r"(\s*NBDirsMax\s*\]\s*\)\s*)\*\s*NBDirsMax\s*\)",
            r"\1)",
            text,
        )
        # Fix broken malloc line: sizeof(double [ NBDirsMax])*NBDirsMax) -> sizeof(double [NBDirsMax]))
        text = re.sub(
            r"sizeof\s*\(\s*double\s*\[\s*NBDirsMax\s*\]\s*\)\s*\*\s*NBDirsMax\s*\)",
            "sizeof(double [NBDirsMax]))",
            text,
        )
        # Zeroing: xb[ii1][nd] -> xb[nd][ii1]
        text = re.sub(r"\bxb\[ii1\]\[nd\]\s*=\s*0\.0", "xb[nd][ii1] = 0.0", text)
        text = re.sub(r"\byb\[ii1\]\[nd\]\s*=\s*0\.0", "yb[nd][ii1] = 0.0", text)
        # xxb[nd] = (xxb+i)[nd] -> memcpy
        text = re.sub(
            r"xxb\[nd\]\s*=\s*\(xxb\+i\)\[nd\]",
            "memcpy(xxb[nd], (xxb+i)[nd], NBDirsMax*sizeof(double))",
            text,
        )
        text = re.sub(
            r"xxb\[nd\]\s*=\s*\((xxb\+[^)]+)\)\[nd\]",
            r"memcpy(xxb[nd], (\1)[nd], NBDirsMax*sizeof(double))",
            text,
        )
        # xb[nd] = (xb+...)[nd]
        text = re.sub(r"xb\[nd\]\s*=\s*\((xb\+[^)]+)\)\[nd\]", r"xb[nd] = (\1)[nd]", text)
        text = re.sub(r"xb\[nd\]\s*=\s*\((xb\+tincx)\)\[nd\]", r"xb[nd] = (xb+tincx)[nd]", text)
        # (*yb)[nd] etc. -> ((double *)yb)[nd]
        text = re.sub(r"\(\*yb\)\[nd\]", "((double *)yb)[nd]", text)
        text = re.sub(r"\(\*xb\)\[nd\]", "((double *)xb)[nd]", text)
        text = re.sub(r"\(\*xxb\)\[nd\]", "((double *)xxb)[nd]", text)
        text = re.sub(r"\(\*yyb\)\[nd\]", "((double *)yyb)[nd]", text)
        # xb[1][nd] / yyb[1][nd] in pop block
        text = re.sub(r"\bxb\[1\]\[nd\]", "((double *)xb)[nd]", text)
        text = re.sub(r"\byyb\[1\]\[nd\]", "((double *)yyb)[nd]", text)
        # yb[nd]++, --yb[nd]
        text = re.sub(r"\byb\[nd\]\+\+", "yb[nd]++", text)
        text = re.sub(r"--yb\[nd\]", "--yb[nd]", text)
        if "memcpy" in text and "#include <stdlib.h>" in text and "#include <string.h>" not in text:
            text = text.replace("#include <stdlib.h>", "#include <stdlib.h>\n#include <string.h>", 1)

    if text != original:
        path.write_text(text)
        return True
    return False


def fix_real_bv_array_type_assignment(path: Path) -> bool:
    """
    Fix "assignment to expression with array type" in real _bv.c files.
    When a matrix adjoint is declared as (*Param)[NBDirsMax], Param[nd] has array type
    and Param[nd] = 0.0 is invalid. Replace with a loop that zeros each element.
    Only apply for params that are actually declared as pointer-to-array in this file.
    Revert any previous mistaken replacement when param is declared as flat pointer (double *Ab).
    Returns True if file was modified.
    """
    text = path.read_text()
    original = text
    for param in ARRAY_TYPE_ADJOINT_PARAMS:
        has_pointer_to_array = bool(
            re.search(rf"\(\s*\*\s*{re.escape(param)}\s*\)\s*\[\s*NBDirsMax\s*\]", text)
        )
        # Flat pointer declaration: "double *Ab" or "float *Ab" (param as separate token)
        has_flat_pointer = bool(
            re.search(rf"(?:double|float)\s+\*\s*{re.escape(param)}\s*[,)]", text)
        )
        if has_flat_pointer:
            # Revert mistaken loop back to simple assignment
            pattern = (
                rf"(\s*)\{{\s*int _ii;\s*for \(_ii = 0;\s*_ii < NBDirsMax;\s*_ii\+\+\)\s*"
                + re.escape(param)
                + r"\[nd\]\[_ii\] = 0\.0;\s*\}"
            )
            text = re.sub(pattern, rf"\1{param}[nd] = 0.0;", text, flags=re.MULTILINE)
        elif has_pointer_to_array:
            # Replace Param[nd] = 0.0 with loop
            pattern = rf"^(\s*){re.escape(param)}\[nd\]\s*=\s*0\.0\s*;"
            replacement = (
                r"\1{ int _ii; for (_ii = 0; _ii < NBDirsMax; _ii++) "
                + param
                + r"[nd][_ii] = 0.0; }"
            )
            text = re.sub(pattern, replacement, text, flags=re.MULTILINE)
    if text != original:
        path.write_text(text)
        return True
    return False


def fix_complex_bv_void_casts_in_dir(src_dir: Path) -> list[str]:
    """
    Apply void* cast fixes to all cblas_c*_bv.c and cblas_z*_bv.c in src_dir.
    Returns the list of filenames that were modified.
    """
    modified: list[str] = []
    if not src_dir.is_dir():
        return modified
    for path in sorted(src_dir.glob("cblas_c*_bv.c")):
        if fix_complex_bv_file(path, is_single=True):
            modified.append(path.name)
    for path in sorted(src_dir.glob("cblas_z*_bv.c")):
        if fix_complex_bv_file(path, is_single=False):
            modified.append(path.name)
    return modified


def fix_real_bv_array_type_in_dir(src_dir: Path) -> list[str]:
    """
    Apply array-type assignment fix to all real _bv.c (cblas_d*_bv.c, cblas_s*_bv.c) in src_dir.
    Returns the list of filenames that were modified.
    """
    modified: list[str] = []
    if not src_dir.is_dir():
        return modified
    for path in sorted(src_dir.glob("cblas_d*_bv.c")) + sorted(src_dir.glob("cblas_s*_bv.c")):
        if fix_real_bv_array_type_assignment(path):
            modified.append(path.name)
    return modified


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    if len(sys.argv) > 1:
        src_dir = Path(sys.argv[1])
        if not src_dir.is_absolute():
            src_dir = (script_dir / src_dir).resolve()
    else:
        src_dir = (script_dir / "out_cblas" / "src").resolve()
        if not src_dir.exists():
            src_dir = script_dir

    if not src_dir.is_dir():
        print(f"Error: {src_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    modified_complex = fix_complex_bv_void_casts_in_dir(src_dir)
    modified_real = fix_real_bv_array_type_in_dir(src_dir)
    if modified_complex:
        print("Fixed void* casts in:", ", ".join(modified_complex))
    if modified_real:
        print("Fixed array-type assignment in:", ", ".join(modified_real))
    if not modified_complex and not modified_real:
        print("No _bv.c files needed changes (or none found).")


if __name__ == "__main__":
    main()
